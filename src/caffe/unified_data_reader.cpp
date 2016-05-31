#include <boost/thread.hpp>
#include <map>
#include <string>
#include <vector>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/convert.hpp>
#include <boost/convert/lexical_cast.hpp>

#include "caffe/unified_data_reader.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/util/rng.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

using boost::weak_ptr;

map<const string, weak_ptr<UnifiedDataReader::Body> > UnifiedDataReader::bodies_;
static boost::mutex unified_bodies_mutex_;

UnifiedDataReader::UnifiedDataReader(const LayerParameter& param) : queue_pair_(new QueuePair(param.unified_data_param().batch_size())) {
  // Get or create a body to share when training
  if (param.phase() == TRAIN) {
    boost::mutex::scoped_lock lock(unified_bodies_mutex_);
    string key = source_key(param);
    weak_ptr<Body>& weak = bodies_[key];
    body_ = weak.lock();
    if (!body_) {
      if (param.unified_data_param().category_base_sampling()) {
        CHECK(param.unified_data_param().has_category_based_sampling_data_param()) << "Should contain category_based_sampling_data_param";
        body_.reset(new Category_Based_Sample_Body(param));
        bodies_[key] = weak_ptr<Body>(body_);
      }
      else {
        body_.reset(new Body(param));
        bodies_[key] = weak_ptr<Body>(body_);
      }
    }
  }
  // No share on the testing
  else {
    if (param.unified_data_param().category_base_sampling()) {
      CHECK(param.unified_data_param().has_category_based_sampling_data_param()) << "Should contain category_based_sampling_data_param";
      body_.reset(new Category_Based_Sample_Body(param));
    }
    else {
      body_.reset(new Body(param));
    }
  }
  body_->new_queue_pairs_.push(queue_pair_);
}

UnifiedDataReader::~UnifiedDataReader() {
  body_.reset();
  if (body_->param_.phase() == TRAIN) {
    string key = source_key(body_->param_);
    boost::mutex::scoped_lock lock(unified_bodies_mutex_);
    if (bodies_[key].expired()) {
      bodies_.erase(key);
    }
  }
}

void UnifiedDataReader::Body::ShuffleSamples() {
  caffe::rng_t* prefetch_rng =
    static_cast<caffe::rng_t*>(prefetch_rng_->generator());
  shuffle(samples_.begin(), samples_.end(), prefetch_rng);
}

UnifiedDataReader::QueuePair::QueuePair(int size) {
  // Initialize the free queue with requested number of datums
  for (int i = 0; i < size; ++i) {
    UnifiedSample* psample = new UnifiedSample();
    psample->image_data = new vector<ImageData>();
    free_.push(psample);
  }
}

UnifiedDataReader::QueuePair::~QueuePair() {
  UnifiedSample* psample;
  while (free_.try_pop(&psample)) {
    delete psample->image_data;
    delete psample;
  }
  while (full_.try_pop(&psample)) {
    delete psample->image_data;
    delete psample;
  }
}

//
UnifiedDataReader::Body::Body(const LayerParameter& param)
    : param_(param),
      new_queue_pairs_() {
  // check load in memory and shuffle
  load_in_memory_ = param.unified_data_param().load_in_memory();
  shuffle_ = param.unified_data_param().shuffle();
  source_ = param_.unified_data_param().source();
  if (shuffle_) CHECK(load_in_memory_) << "load_in_memory must be set to true if you need shuffling.";
  if (load_in_memory_) {
    CHECK(samples_.empty()) << "samples are not empty -- read files should happen for only once";
    string line;
    size_t num = 0;
    sample_id_ = 0;

    LOG(INFO) << "Opening file " << source_;
    data_stream_.open(source_.c_str());

    while (getline(data_stream_, line)) {
      if (line == "") break;
      ++num;
      samples_.resize(num);
      auto& p_sample = samples_[num - 1];
      p_sample.reset(new UnifiedSample);
      auto& sample = *p_sample;
      vector <string> ps;
      boost::split(ps, line, boost::is_any_of("\t"));
      try {
        sample.label = boost::lexical_cast<float>(ps[ps.size() - 1]);
      }
      catch (boost::bad_lexical_cast const&) {
        CHECK(false) << "the label should be float, part: " << ps[ps.size() - 1] << ", line: " << line;
      }
      sample.data.resize(ps.size() - 1);
      for (size_t i = 0; i < ps.size() - 1; ++i){
        sample.data[i] = ps[i];
      }
    }

    if (shuffle_) {
      const unsigned int prefetch_rng_seed = caffe_rng_rand();
      prefetch_rng_.reset(new Caffe::RNG(prefetch_rng_seed));
      ShuffleSamples();
    }
    CHECK(!samples_.empty()) << "no samples were read from input file: " << source_;
    LOG(INFO) << "A total of " << samples_.size() << " samples.";
  }
  if (param.unified_data_param().category_base_sampling() == false)
    StartInternalThread();
}

UnifiedDataReader::Body::~Body() {
  StopInternalThread();
}

void UnifiedDataReader::Body::InternalThreadEntry() {
  vector<shared_ptr<QueuePair> > qps;
  try {
    int solver_count = param_.phase() == TRAIN ? Caffe::solver_count() : 1;
    // To ensure deterministic runs, only start running once all solvers
    // are ready. But solvers need to peek on one item during initialization,
    // so read one item, then wait for the next solver.
    for (int i = 0; i < solver_count; ++i) {
      shared_ptr<QueuePair> qp(new_queue_pairs_.pop());
      read_one(qp.get());
      qps.push_back(qp);
    }
    // Main loop
    while (!must_stop()) {
      for (int i = 0; i < solver_count; ++i) {
        read_one(qps[i].get());
      }
      // Check no additional readers have been created. This can happen if
      // more than one net is trained at a time per process, whether single
      // or multi solver. It might also happen if two data layers have same
      // name and same source.
      CHECK_EQ(new_queue_pairs_.size(), 0);
    }
  }
  catch (boost::thread_interrupted&) {
    // Interrupted exception is expected on shutdown
  }
}

void UnifiedDataReader::Body::read_one(QueuePair* qp) {
  UnifiedSample* pstr = qp->free_.pop();
  if (load_in_memory_)
    PrefetchOneInMemory(pstr);
  else
    PrefetchOneFromReadingFile(pstr);
  PreloadIUB(pstr);
  qp->full_.push(pstr);
}

void UnifiedDataReader::Body::PrefetchOneInMemory(UnifiedSample* pstr){
  if (sample_id_ == samples_.size()){
    if (shuffle_) {
      ShuffleSamples();
    }
    sample_id_ = 0;
  }
  pstr->data.resize(samples_[sample_id_]->data.size());
  for (size_t i = 0; i < samples_[sample_id_]->data.size(); ++i){
    pstr->data[i] = samples_[sample_id_]->data[i];
  }
  pstr->label = (float)samples_[sample_id_]->label;
  sample_id_++;
}

void UnifiedDataReader::Body::PrefetchOneFromReadingFile(UnifiedSample* pstr){
  string line;
  if (!std::getline(data_stream_, line)){
    data_stream_.close();
    data_stream_.open(source_.c_str());
    CHECK(std::getline(data_stream_, line)) << "can't read file";
  }
  vector <string> ps;
  boost::split(ps, line, boost::is_any_of("\t"));
  try {
    pstr->label = boost::lexical_cast<float>(ps[ps.size() - 1]);
  }
  catch (boost::bad_lexical_cast const&) {
    CHECK(false) << "the label should be float, part: " << ps[ps.size() - 1] << ", line: " << line;
  }
  pstr->data.resize(ps.size() - 1);
  for (size_t i = 0; i < ps.size() - 1; ++i){
    pstr->data[i] = ps[i];
  }
}

void UnifiedDataReader::Body::PreloadIUB(UnifiedSample* pstr) {
  int data_size = pstr->data.size();
  // Create ImageDBs if not created
  if (imageDBs.size() != data_size) {
    imageDBs.resize(data_size);
    for (int i = 0; i < data_size; i++) {
      if (param_.unified_data_param().data_desc(i).type() == DataDescriptionParameter::IMAGEIUB) {
        imageDBs[i] = IImageDB::CreateInstanceAndInit(param_.unified_data_param().data_desc(i));
      }
    }
  }

  // Preload images for all targets that are IUB files
  //pstr->image_data = new vector<ImageData>(data_size);
  (pstr->image_data)->resize(data_size);
  for (int i = 0; i < data_size; i++) {
    if (param_.unified_data_param().data_desc(i).type() == DataDescriptionParameter::IMAGEIUB) {
      imageDBs[i]->Read(pstr->data[i], (*(pstr->image_data))[i]);
    }
  }
}

UnifiedDataReader::Category_Based_Sample_Body::Category_Based_Sample_Body(const LayerParameter& param) : Body(param) {
  samples_per_class_ = param.unified_data_param().category_based_sampling_data_param().samples_per_class();
  CHECK(samples_per_class_ > 0) << "at least sample one from each category.";
  int solver_count = param_.phase() == TRAIN ? Caffe::solver_count() : 1;
  batch_size_ = param.unified_data_param().batch_size() * solver_count;
  if (param.unified_data_param().shuffle()) {
    const unsigned int prefetch_rng_seed = caffe_rng_rand();
    prefetch_rng_.reset(new Caffe::RNG(prefetch_rng_seed));
  }
  const unsigned int sample_rng_seed = caffe_rng_rand();
  sample_rng_.reset(new Caffe::RNG(sample_rng_seed));
  class_sample_list_.clear();
  sample_idx_in_ori_list_.clear();
  idx_to_class_.clear();
  int label_idx = 0;
  for (size_t i = 0; i < samples_.size(); ++i)
  {
    auto it = class_sample_list_.find(samples_[i]->label);
    if (it != class_sample_list_.end()) it->second.push_back(int(i));
    else class_sample_list_.insert(std::make_pair(samples_[i]->label, std::vector<int>(1, int(i))));
  }
  for (auto it = class_sample_list_.begin(); it != class_sample_list_.end();)
  {
    if (samples_per_class_ > 0 && it->second.size() < (size_t)samples_per_class_)
    {
      erased_class_.push_back(it->first);
      it = class_sample_list_.erase(it);
    }
    else
    {
      for (int i = 0; i < it->second.size(); ++i)
      {
        int id = sample_idx_in_ori_list_.size();
        sample_idx_in_ori_list_.push_back(it->second[i]);
        it->second[i] = id;
      }
      idx_to_class_[label_idx++] = it->first;
      class_sample_list_len_.push_back(it->second.size());
      CHECK(label_idx == class_sample_list_len_.size());
      ++it;
    }
  }

  // LOG how many classies are removed
  LOG(INFO) << erased_class_.size() << " out of " << erased_class_.size() + class_sample_list_.size() << " classes are filtered out for they have less samples than " << samples_per_class_;

  if (batch_size_ >(int)sample_idx_in_ori_list_.size())
  {
    LOG(FATAL) << "Not enough number of data for the batch size";
  }
  if (samples_per_class_ > 0 && batch_size_ > (int)class_sample_list_.size() * samples_per_class_)
  {
    LOG(FATAL) << "Not enough number of classes for the batch size";
  }
  StartInternalThread();
}

void UnifiedDataReader::Category_Based_Sample_Body::InternalThreadEntry() {
  vector<shared_ptr<QueuePair> > qps;
  try {
    int solver_count = param_.phase() == TRAIN ? Caffe::solver_count() : 1;
    // To ensure deterministic runs, only start running once all solvers
    // are ready. But solvers need to peek on one item during initialization,
    // so read one item, then wait for the next solver.

    for (int i = 0; i < solver_count; ++i) {
      shared_ptr<QueuePair> qp(new_queue_pairs_.pop());      
      
      // you need to read one to let layer setup thread to continue working, otherwise, there's deadlock
      QueuePair* pqp = qp.get();
      UnifiedSample* pstr = pqp->free_.pop();
      auto& vec = class_sample_list_[idx_to_class_[0]];
      int sampleid = vec[0];
      pstr->data.resize(samples_[sample_idx_in_ori_list_[sampleid]]->data.size());
      for (size_t j = 0; j < pstr->data.size(); ++j){
        pstr->data[j] = samples_[sample_idx_in_ori_list_[sampleid]]->data[j];
      }
      pstr->label = samples_[sample_idx_in_ori_list_[sampleid]]->label;
      pqp->full_.push(pstr);

      qps.push_back(qp);
    }

    /*
    // Empty the queues for we sample exactly one batch_size samples
    for (int i = 0; i < solver_count; ++i) {
      QueuePair* qp = qps[i].get();
      UnifiedSample* pstr = qp->full_.pop();
      CHECK_EQ(qp->full_.size(), 0) << "Should be empty queue.";
      qp->free_.push(pstr);
    }
    */

    // Main loop
    while (!must_stop()) {
      caffe::rng_t* sample_rng = static_cast<caffe::rng_t*>(sample_rng_->generator());
      RandomPermutation class_perm(class_sample_list_len_);
      for (int num = 0; num < batch_size_;) {
        // first random get next class
        float classid = idx_to_class_[class_perm.GetNext(sample_rng)];
        auto& vec = class_sample_list_[classid];
        // second permulate with the samples of a class
        RandomPermutation sample_in_class_perm((int)(vec.size() - 1));
        for (int numc = 0; numc < samples_per_class_ && num < batch_size_; ++numc){
          int solver_idx = num % solver_count;
          QueuePair* qp = qps[solver_idx].get();
          UnifiedSample* pstr = qp->free_.pop();
          int sampleid = vec[sample_in_class_perm.GetNext(sample_rng)];
          pstr->data.resize(samples_[sample_idx_in_ori_list_[sampleid]]->data.size());
          for (size_t i = 0; i < pstr->data.size(); ++i){
            pstr->data[i] = samples_[sample_idx_in_ori_list_[sampleid]]->data[i];
          }
          pstr->label = samples_[sample_idx_in_ori_list_[sampleid]]->label;
          qp->full_.push(pstr);
          num++;
        }
      }
      // Check no additional readers have been created. This can happen if
      // more than one net is trained at a time per process, whether single
      // or multi solver. It might also happen if two data layers have same
      // name and same source.
      CHECK_EQ(new_queue_pairs_.size(), 0);
    }
  }
  catch (boost::thread_interrupted&) {
    // Interrupted exception is expected on shutdown
  }
}


}  // namespace caffe
