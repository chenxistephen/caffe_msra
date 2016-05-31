#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <ppl.h>

#include "caffe/layer.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/util/rng.hpp"
#include "caffe/layers/jointsimilarity_layer.hpp"
#include "Base.h"
#include "FastOperations.h"

namespace caffe {
    template <typename Dtype>
    void JointSimilarityLayer<Dtype>::LayerSetUp(
        const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
        Layer<Dtype>::LayerSetUp(bottom, top);
        pairwise_ = this->layer_param_.joint_similarity_param().pairwise();
        column_sample_num_ = this->layer_param_.joint_similarity_param().column_sample_num();
        eps_ = this->layer_param_.joint_similarity_param().eps();
        if (pairwise_) {
            CHECK(bottom.size() == 2) << "For pairwise similarity, there should be 2 inputs.";
            CHECK_EQ(bottom[1]->count(), bottom[0]->count()) << "The 2 input should have the same dim.";
        }
        CHECK(column_sample_num_ <= bottom[0]->num()) << "We can't sample more columns than batchsize.";
        CHECK(column_sample_num_ > 0) << "At least sample one column";
        column_sample_rng_.reset(new Caffe::RNG(caffe_rng_rand()));
    }

    template <typename Dtype>
    void JointSimilarityLayer<Dtype>::Reshape(
        const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
        // Deducing top shape of output
        vector<int> top_shape;
        top_shape.push_back(bottom[0]->num());
        top_shape.push_back(column_sample_num_);
        top[0]->Reshape(top_shape);
        // For sample_mat, sample_mat[i][0] indicates how many columns are sampled for i
        vector<int> sample_mat_shape;
        sample_mat_shape.push_back(bottom[0]->num());
        sample_mat_shape.push_back(bottom[0]->num());
        sample_mat0_.Reshape(sample_mat_shape);
        sample_mat1_.Reshape(sample_mat_shape);
        // prepare space for derivative for each pair
        vector<int> der_shape;
        der_shape.push_back(bottom[0]->count() * column_sample_num_);
        derivative_mat0_.Reshape(der_shape);
        derivative_mat1_.Reshape(der_shape);
    }

    template <typename Dtype>
    Dtype L2JointSimilarityLayer<Dtype>::ComputeSimilarity_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, int dim, const Dtype eps) {
        Dtype res = 0;
        for (int i = 0; i < dim; i++)
            res += (bottom_row0[i] - bottom_row1[i]) * (bottom_row0[i] - bottom_row1[i]);
        return -res;
    }

    template <typename Dtype>
    void L2JointSimilarityLayer<Dtype>::PropagateDown_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, Dtype top_diff, Dtype* bottom_diff_row, int dim, const Dtype eps) {
        for (int i = 0; i < dim; i++)
            bottom_diff_row[i] += (2 * top_diff * (bottom_row1[i] - bottom_row0[i]));
    }

    template <typename Dtype>
    Dtype CosineJointSimilarityLayer<Dtype>::ComputeSimilarity_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, int dim, const Dtype eps) {
        Dtype dotproduct_data = caffe_cpu_dot(dim, bottom_row0, bottom_row1);
        Dtype norm_data0 = sqrt(caffe_cpu_dot(dim, bottom_row0, bottom_row0));
        Dtype norm_data1 = sqrt(caffe_cpu_dot(dim, bottom_row1, bottom_row1));
        return 0.5 * dotproduct_data / (norm_data0 * norm_data1 + eps) + 0.5;
    }

    template <typename Dtype>
    void CosineJointSimilarityLayer<Dtype>::PropagateDown_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, Dtype top_diff, Dtype* bottom_diff_row, int dim, const Dtype eps) {
        Dtype dotproduct_data = caffe_cpu_dot(dim, bottom_row0, bottom_row1);
        Dtype norm_data0 = sqrt(caffe_cpu_dot(dim, bottom_row0, bottom_row0));
        Dtype norm_data1 = sqrt(caffe_cpu_dot(dim, bottom_row1, bottom_row1));
        Dtype alpha = 0.5 / std::max(norm_data0 * norm_data1, eps_);
        Dtype beta = -0.5 * dotproduct_data / std::max(pow(norm_data0, Dtype(3.0)) * norm_data1, eps_);
        for (int i = 0; i < dim; i++) {
            bottom_diff_row[i] += (top_diff * (alpha * bottom_row1[i] + beta * bottom_row0[i]));
        }
    }

    template <typename Dtype>
    void JointSimilarityLayer<Dtype>::Sample_pairs(int num, int dim) {
        Dtype* sample_mat_array0 = sample_mat0_.mutable_cpu_data();
        Dtype* sample_mat_array1 = sample_mat1_.mutable_cpu_data();
        caffe_memset(num * num * sizeof(Dtype), 0, sample_mat_array0);
        caffe_memset(num * num * sizeof(Dtype), 0, sample_mat_array1);

        // Only compute matched pair
        if (column_sample_num_ == 1) {
            for (int i = 0; i < num; i++) {
                sample_mat_array0[i * num] = i + 1;
                sample_mat_array1[i * num] = i + 1;
            }
        }
        // Compute every n * n pairs
        else if (column_sample_num_ == num) {
            for (int i = 0; i < num; i++) {
                // place the postive sample at the first column
                sample_mat_array0[i * num] = i + 1;
                sample_mat_array1[i * num + i] = i * num + 1;
                for (int j = 1; j < num; j++) {
                    if (j <= i) {
                        sample_mat_array0[i * num + j] = j + 1 - 1;
                        sample_mat_array1[(j - 1) * num + i] = i * num + j + 1;
                    }
                    else{
                        sample_mat_array0[i * num + j] = j + 1;
                        sample_mat_array1[j * num + i] = i * num + j + 1;
                    }
                }
            }
        }
        // Sample k pairs, dist[i, i] must be sampled
        else {
            for (int i = 0; i < num; i++) {
                // ith bottom0 and ith bottom1 must be sampled for it is considered as positive pair
                sample_mat_array0[i * num] = i + 1;
                for (int n = 0; n < num; n++) {
                    if (sample_mat_array1[i * num + n] == 0) {
                        sample_mat_array1[i * num + n] = i * column_sample_num_ + 1;
                        break;
                    }
                }
                RandomPermutation rand_perm(num - 2);
                for (int k = 1; k < column_sample_num_; k++) {
                    //int neg_idx = rand_perm.GetNext(col_rng);
                    int neg_idx = k - 1;
                    neg_idx = (neg_idx < i) ? neg_idx : neg_idx + 1;
                    sample_mat_array0[i * num + k] = neg_idx + 1;
                    for (int n = 0; n < num; n++) {
                        if (sample_mat_array1[neg_idx * num + n] == 0) {
                            sample_mat_array1[neg_idx * num + n] = i * column_sample_num_ + k + 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    template <typename Dtype>
    void JointSimilarityLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
        const Dtype* bottom_data0 = bottom[0]->cpu_data();
        const Dtype* bottom_data1 = pairwise_ ? bottom[1]->cpu_data() : bottom[0]->cpu_data();
        Dtype* dist = top[0]->mutable_cpu_data();
        Dtype* sample_mat_array0 = sample_mat0_.mutable_cpu_data();
        Dtype* sample_mat_array1 = sample_mat1_.mutable_cpu_data();
        int num = bottom[0]->num();
        int dim = bottom[0]->count() / bottom[0]->num();
        caffe_memset(num * column_sample_num_ * sizeof(Dtype), 0, dist);
        Sample_pairs(num, dim);
        concurrency::parallel_for(int(0), num, [&](int t) {
            const Dtype* bottom_row0 = bottom_data0 + t * dim;
            for (int k = 0; k < column_sample_num_; k++) {
                const Dtype* bottom_row1 = bottom_data1 + int(sample_mat_array0[t * num + k] - 1) * dim;
                dist[t * column_sample_num_ + k] = ComputeSimilarity_cpu(bottom_row0, bottom_row1, dim, eps_);
            }
        });    
    }

    template <typename Dtype>
    void JointSimilarityLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
        const Dtype* bottom_data0 = bottom[0]->cpu_data();
        const Dtype* bottom_data1 = pairwise_ ? bottom[1]->cpu_data() : bottom[0]->cpu_data();
        const Dtype* sample_mat_array0 = sample_mat0_.mutable_cpu_data();
        const Dtype* sample_mat_array1 = sample_mat1_.mutable_cpu_data();
        const Dtype* top_diff = top[0]->cpu_diff();
        Dtype* bottom_diff0 = propagate_down[0] ? bottom[0]->mutable_cpu_diff() : nullptr;
        Dtype* bottom_diff1 = propagate_down[1] ? (pairwise_ ? bottom[1]->mutable_cpu_diff() : bottom[0]->mutable_cpu_diff()) : nullptr;
        if (propagate_down[0])
            caffe_set(bottom[0]->count(), Dtype(0), bottom_diff0);
        if (propagate_down[1] && pairwise_)
            caffe_set(bottom[1]->count(), Dtype(0), bottom_diff1);
        int num = bottom[0]->num();
        int dim = bottom[0]->count() / bottom[0]->num();
        // Only compute matched pair
        // [TODO] haven't checked !pairwise correctness
        if (propagate_down[0]) {
            concurrency::parallel_for(int(0), num, [&](int k0) {
                const Dtype* bottom_row_target = bottom_data0 + k0 * dim;
                Dtype* bottom_diff_row = bottom_diff0 + k0 * dim;
                for (int k1 = 0; k1 < column_sample_num_; ++k1) {
                    const Dtype* bottom_row_sibling = bottom_data1 + int(sample_mat_array0[k0 * num + k1] - 1) * dim;
                    PropagateDown_cpu(bottom_row_target, bottom_row_sibling, top_diff[k0 * column_sample_num_ + k1], bottom_diff_row, dim, eps_);
                }
            });
        }
        if (propagate_down[1]) {
            concurrency::parallel_for(int(0), num, [&](int k0) {
                const Dtype* bottom_row_target = bottom_data1 + k0 * dim;
                Dtype* bottom_diff_row = bottom_diff1 + k0 * dim;
                for (int k1 = 0; k1 < num && sample_mat_array1[k0 * num + k1] > 0; ++k1) {
                    int siblingidx = int(sample_mat_array1[k0 * num + k1] - 1) / column_sample_num_;
                    const Dtype* bottom_row_sibling = bottom_data0 + siblingidx * dim;
                    PropagateDown_cpu(bottom_row_target, bottom_row_sibling, top_diff[int(sample_mat_array1[k0 * num + k1] - 1)], bottom_diff_row, dim, eps_);
                }
            });
        }
    }



#ifdef CPU_ONLY
    STUB_GPU(L2JointSimilarityLayer);
    STUB_GPU(CosineJointSimilarityLayer);
#endif

INSTANTIATE_CLASS(L2JointSimilarityLayer);
REGISTER_LAYER_CLASS(L2JointSimilarity);
INSTANTIATE_CLASS(CosineJointSimilarityLayer);
REGISTER_LAYER_CLASS(CosineJointSimilarity);
}  // namespace caffe