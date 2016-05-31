#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <ppl.h>

#include "caffe/layer.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/layers/nca_loss_layer.hpp"
#include "caffe/util/FastOperations.h"

namespace caffe {

    template <typename Dtype>
    void NcaLossLayer<Dtype>::Reshape(
        const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
        LossLayer<Dtype>::Reshape(bottom, top);
        pairwise_ = this->layer_param_.nca_loss_param().pairwise();
        if (pairwise_ == true) {
            CHECK_EQ(bottom[1]->channels(), bottom[0]->channels());
            CHECK_EQ(bottom[1]->height(), bottom[0]->height());
            CHECK_EQ(bottom[1]->width(), bottom[0]->width());
        }
        else {
            CHECK_EQ(bottom[1]->count(2), 1) << "input dimension is incorrect, this should be the label, where only have one dim";
        }
    }

    template <typename Dtype>
    void NcaLossLayer<Dtype>::CalculateDist_cpu(const Dtype* bottom_data, const Dtype* bottom_data1, int num, int dim, std::vector<Dtype>& dist) {
        dist.assign(num * num, Dtype(0));
        if (pairwise_) {
            // for pairwise, dist[i, j] is dist(left_i, right_j), dist[j, i] is dist(left_j, right_i), they are different
            int num2 = (num + 1) / 2;
            concurrency::parallel_for(int(0), num2 * num2, [&](int t) {
                int b0 = (t / num2), b1 = (t % num2);
                b0 *= 2; b1 *= 2;
                int e0 = std::min(b0 + 2, num), e1 = std::min(b1 + 2, num);
                for (int k0 = b0; k0 < e0; ++k0) {
                    const Dtype* bottom_row0 = bottom_data + k0 * dim;
                    for (int k1 = b1; k1 < e1; ++k1) {
                        const Dtype* bottom_row1 = bottom_data1 + k1 * dim;
                        Dtype sum = ns_base::ns_sse::EuclideanDistA(bottom_row0, bottom_row1, dim);
                        dist[k0 * num + k1] = exp(-sum);
                    }
                }
            });
        }
        else {
            // dist[i, j] == dist[j, i], so only half of the matrix need computing
            int num2 = (num + 1) / 2;
            int num22 = (num2 + 1) / 2;
            concurrency::parallel_for(int(0), (num2 + 1) * num22, [&](int t) {
                int b0 = (t / num22), b1 = (t % num22);
                if (b0 > b1) {
                    if (2 * b1 == num2 - 1) return;
                    b0 = num2 - b0;
                    b1 = num2 - 1 - b1;
                }
                b0 *= 2; b1 *= 2;
                int e0 = std::min(b0 + 2, num), e1 = std::min(b1 + 2, num);
                for (int k0 = b0; k0 < e0; ++k0) {
                    const Dtype* bottom_row0 = bottom_data + k0 * dim;
                    for (int k1 = std::max(b1, k0 + 1); k1 < e1; ++k1) {
                        const Dtype* bottom_row1 = bottom_data + k1 * dim;
                        Dtype sum = ns_base::ns_sse::EuclideanDistA(bottom_row0, bottom_row1, dim);
                        dist[k1 * num + k0] = dist[k0 * num + k1] = exp(-sum);
                    }
                }
            });
        }
    }

    template <typename Dtype>
    Dtype NcaLossLayer<Dtype>::ForwardInternal(const Dtype* dist, const Dtype* bottom_label, int num) {
        int sum_num = 0;
        Dtype sum = 0;
        for (int k = 0; k < num; ++k) {
            bool flag = false;
            Dtype pos = small_const_, neg = Dtype(1e6) * small_const_;
            for (int i = 0; i < num; ++i) {
                if (pairwise_ || i != k) {
                    if (i == k || (!pairwise_ && bottom_label[i] == bottom_label[k])) {
                        flag = true;
                        pos += dist[k * num + i];
                    }
                    else {
                        neg += dist[k * num + i];
                    }
                }
            }
            if (flag) {
                sum_num += 1;
                sum += -log(pos / (pos + neg));
            }
            if (pairwise_)
            {
                pos = small_const_, neg = Dtype(1e6) * small_const_;
                for (int i = 0; i < num; ++i) {
                    if (i == k) {
                        pos += dist[i * num + k];
                    }
                    else {
                        neg += dist[i * num + k];
                    }
                }
                sum_num += 1;
                sum += -log(pos / (pos + neg));
            }
        }
        return sum / sum_num;
    }

    template <typename Dtype>
    void NcaLossLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
        const Dtype* bottom_data = bottom[0]->cpu_data();
        const Dtype* bottom_data1 = pairwise_ ? bottom[1]->cpu_data() : nullptr;
        const Dtype* bottom_label = !pairwise_ ? bottom[1]->cpu_data() : nullptr;
        int num = bottom[0]->num();
        int dim = bottom[0]->count() / bottom[0]->num();
        std::vector<Dtype> dist;
        CalculateDist_cpu(bottom_data, bottom_data1, num, dim, dist);
        top[0]->mutable_cpu_data()[0] = ForwardInternal(dist.data(), bottom_label, num);
    }

    template <typename Dtype>
    void NcaLossLayer<Dtype>::BackwardInternal(const Dtype* dist, const Dtype* bottom_label, int num, Dtype mult, std::vector<Dtype>& dist_grad) {
        int sum_num = 0;
        dist_grad.assign(num * num, Dtype(0));
        for (int k = 0; k < num; ++k) {
            bool flag = false;
            Dtype pos = small_const_, neg = Dtype(1e6) * small_const_;
            for (int i = 0; i < num; ++i) {
                if (pairwise_ || i != k) {
                    if (i == k || (!pairwise_ && bottom_label[i] == bottom_label[k])) {
                        flag = true;
                        pos += dist[k * num + i];
                    }
                    else {
                        neg += dist[k * num + i];
                    }
                }
            }
            if (flag) {
                sum_num += 1;
                Dtype pos_mult = mult * (1.0 / (pos + neg) - 1.0 / pos);
                Dtype neg_mult = mult / (pos + neg);
                for (int i = 0; i < num; ++i) {
                    if (pairwise_ || i != k) {
                        if (i == k || (!pairwise_ && bottom_label[i] == bottom_label[k])) {
                            dist_grad[k * num + i] += pos_mult;
                        }
                        else {
                            dist_grad[k * num + i] += neg_mult;
                        }
                    }
                }
            }
            if (pairwise_)
            {
                pos = small_const_, neg = Dtype(1e6) * small_const_;
                for (int i = 0; i < num; ++i) {
                    if (i == k) {
                        flag = true;
                        pos += dist[i * num + k];
                    }
                    else {
                        neg += dist[i * num + k];
                    }
                }
                sum_num += 1;
                Dtype pos_mult = mult * (1.0 / (pos + neg) - 1.0 / pos);
                Dtype neg_mult = mult / (pos + neg);
                for (int i = 0; i < num; ++i) {
                    if (pairwise_ || i != k) {
                        if (i == k) {
                            dist_grad[i * num + k] += pos_mult;
                        }
                        else {
                            dist_grad[i * num + k] += neg_mult;
                        }
                    }
                }
            }
        }
        if (sum_num > 0) {
            Dtype rev_sum_num = Dtype(1) / Dtype(sum_num);
            for (int i = 0; i < num * num; ++i) dist_grad[i] *= rev_sum_num;
        }
    }

    template <typename Dtype>
    void NcaLossLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
        if (propagate_down[1] && !pairwise_) {
            LOG(FATAL) << this->type()
                << " Layer cannot backpropagate to label inputs.";
        }
        if (propagate_down[0] || propagate_down[1]) {
            const Dtype* bottom_data = bottom[0]->cpu_data();
            const Dtype* bottom_data1 = pairwise_ ? bottom[1]->cpu_data() : nullptr;
            const Dtype* bottom_label = !pairwise_ ? bottom[1]->cpu_data() : nullptr;
            Dtype* bottom_diff = propagate_down[0] ? bottom[0]->mutable_cpu_diff() : nullptr;
            Dtype* bottom_diff1 = propagate_down[1] ? bottom[1]->mutable_cpu_diff() : nullptr;
            int num = bottom[0]->num();
            int dim = bottom[0]->count() / bottom[0]->num();
            caffe_set(bottom[0]->count(), Dtype(0), bottom_diff);
            std::vector<Dtype> dist(num * num, Dtype(0));
            CalculateDist_cpu(bottom_data, bottom_data1, num, dim, dist);
            std::vector<Dtype> dist_grad;
            BackwardInternal(dist.data(), bottom_label, num, top[0]->cpu_diff()[0], dist_grad);
            if (pairwise_) {
                if (propagate_down[0]) {
                    concurrency::parallel_for(int(0), num, [&](int k0) {
                        const Dtype* bottom_row0 = bottom_data + k0 * dim;
                        Dtype* bottom_grad_row0 = bottom_diff + k0 * dim;
                        for (int k1 = 0; k1 < num; ++k1) {
                            const Dtype* bottom_row1 = bottom_data1 + k1 * dim;
                            Dtype pmult = -dist_grad[k0 * num + k1] * dist[k0 * num + k1] * 2;
                            ns_base::ns_sse::MultByConstAndAddU(bottom_row0, pmult, dim, bottom_grad_row0);
                            ns_base::ns_sse::MultByConstAndAddU(bottom_row1, -pmult, dim, bottom_grad_row0);
                        }
                    });
                }
                if (propagate_down[1]) {
                    concurrency::parallel_for(int(0), num, [&](int k1) {
                        const Dtype* bottom_row1 = bottom_data1 + k1 * dim;
                        Dtype* bottom_grad_row1 = bottom_diff1 + k1 * dim;
                        for (int k0 = 0; k0 < num; ++k0) {
                            const Dtype* bottom_row0 = bottom_data + k0 * dim;
                            Dtype pmult = -dist_grad[k0 * num + k1] * dist[k0 * num + k1] * 2;
                            ns_base::ns_sse::MultByConstAndAddU(bottom_row0, -pmult, dim, bottom_grad_row1);
                            ns_base::ns_sse::MultByConstAndAddU(bottom_row1, pmult, dim, bottom_grad_row1);
                        }
                    });
                }
            }
            else {
                concurrency::parallel_for(int(0), num, [&](int k0) {
                    const Dtype* bottom_row0 = bottom_data + k0 * dim;
                    Dtype* bottom_grad_row0 = bottom_diff + k0 * dim;
                    for (int k1 = 0; k1 < num; ++k1) {
                        if (k1 != k0) {
                            const Dtype* bottom_row1 = bottom_data + k1 * dim;
                            Dtype pmult = -(dist_grad[k1 * num + k0] + dist_grad[k0 * num + k1]) * dist[k0 * num + k1] * 2;
                            ns_base::ns_sse::MultByConstAndAddU(bottom_row0, pmult, dim, bottom_grad_row0);
                            ns_base::ns_sse::MultByConstAndAddU(bottom_row1, -pmult, dim, bottom_grad_row0);
                        }
                    }
                });
            }
        }
    }

#ifdef CPU_ONLY
    STUB_GPU(NcaLossLayer);
#endif

INSTANTIATE_CLASS(NcaLossLayer);
REGISTER_LAYER_CLASS(NcaLoss);

}  // namespace caffe
