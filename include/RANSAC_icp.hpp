#pragma once
#include <vector>
#include "opencv2/opencv.hpp"
#include <opencv2/core/core.hpp>
#include <math.h>
#include <stdlib.h>
#include <time.h> 
#include <Eigen/Core>
#include <Eigen/SVD> 
#include <Eigen/Geometry>

// #define DEBUG

class RansacRigidTransform 
{
public:
	RansacRigidTransform(float prob = 0.995, float inlier_p = 0.8, float threshold = 0.02, int m = 6) :_prob(prob), _inlier_p(inlier_p), _threshold(threshold), _m(m)
	{
		srand((unsigned)time(NULL));  // random seed, seed different -> random val different
	};

	~RansacRigidTransform() {};

	int get_N();

	// sample _m points
	void sample_index_lst(const int &min_index, const int &max_index, const int &m, int * const index_lst);

	// calculate rigid transform(ICP)
	void calculate_rigid_transform(const std::vector<Eigen::Vector3d> &src, const std::vector<Eigen::Vector3d> &dst, 
		const int * const index_lst, Eigen::Isometry3d &eu3);

	// get transform distance
	double get_distance(const Eigen::Vector3d &src_pts, const Eigen::Vector3d &dst_pts, const Eigen::Isometry3d &eu3);

	void start_RANSAC(const std::vector<Eigen::Vector3d> &src, const std::vector<Eigen::Vector3d> &dst, 
		std::vector<uchar> &ransac_status, int max_sampling=1000);

public:
	Eigen::Isometry3d res_euc3;

private:
	float _prob, _inlier_p;
	float _threshold;
	int _m;

};



void RansacRigidTransform::sample_index_lst(const int &min_index, const int &max_index, const int &m, int * const index_lst)
{
	// sampling m integer number between [min_index, max_index]
	for (int i = 0; i < m; ++i) {
		int index = rand() % (max_index - min_index + 1) + min_index;
		index_lst[i] = index;
	}
#ifdef DEBUG
	printf("index_lst = \n");
	for (int i = 0; i < m; ++i) {
		printf("%d, ", index_lst[i]);
	}

#endif
};

int RansacRigidTransform::get_N() {
	// calculate sampling times

	int sample_times = static_cast<int>(log(1 - this->_prob) / log(1 - pow(this->_inlier_p, this->_m) + 1e-8));
	return sample_times;

};

void RansacRigidTransform::calculate_rigid_transform(const std::vector<Eigen::Vector3d> &src, const std::vector<Eigen::Vector3d> &dst, 
		const int * const index_lst, Eigen::Isometry3d &eu3)
{
	// ////////////////////////////////////////////////////
	// calculate src-to-dst rigid transform
	// ////////////////////////////////////////////////////

	Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
	Eigen::Vector3d mean_src = Eigen::Vector3d::Zero();
	Eigen::Vector3d mean_dst = Eigen::Vector3d::Zero();
	for(int i = 0; i < _m; ++i)
	{
		int index = index_lst[i];
		mean_src += src[index];
		mean_dst += dst[index];
	}
	mean_src = mean_src / double(_m);
	mean_dst = mean_dst / double(_m);
	for(int i = 0; i < _m; ++i)
	{
		int index = index_lst[i];
		H = H + (src[index] - mean_src) * (dst[index] - mean_dst).transpose();
	}
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV );
	Eigen::Matrix3d R = svd.matrixV() * svd.matrixU().transpose();
	Eigen::Vector3d trans = mean_dst - R * mean_src;
	eu3 = Eigen::Isometry3d::Identity();
	eu3.rotate(R);
	eu3.pretranslate(trans);
	
};


double RansacRigidTransform::get_distance(const Eigen::Vector3d &src_pts, const Eigen::Vector3d &dst_pts, 
const Eigen::Isometry3d &eu3) 
{
	Eigen::Vector3d delta_pos = dst_pts - eu3 * src_pts;
	
	return delta_pos.norm();
};


void RansacRigidTransform::start_RANSAC(const std::vector<Eigen::Vector3d> &src, const std::vector<Eigen::Vector3d> &dst, 
	std::vector<uchar> &ransac_status, int max_sampling)
{
	int num = this->get_N();
	int sampling_num = MIN(num, max_sampling);
	int *index_lst = new int[src.size()];
	int match_size = src.size();
	
	printf("sampling %d times\n", sampling_num);
	printf("matchsize = %d\n", match_size);

	int max_inlier_num = -1;
	for (int i = 0; i < sampling_num; i++)
	{
		sample_index_lst(0, src.size() - 1, _m, index_lst);
		Eigen::Isometry3d euc3;
		this->calculate_rigid_transform(src, dst, index_lst, euc3);
		std::vector<uchar> linshi_ransac_status(match_size);

		// calculate inlier points
		int inlier_num = 0;
		for (int j = 0; j < match_size; ++j) {
			double dist = get_distance(src[j], dst[j], euc3);
			if (dist < _threshold) {
				linshi_ransac_status[j] = 1;
				inlier_num += 1;
			}
			else {
				printf("dist = %f\n", dist);
				linshi_ransac_status[j] = 0;
			}
		}
		if (inlier_num > max_inlier_num) {
			max_inlier_num = inlier_num;
			ransac_status.swap(linshi_ransac_status);
			res_euc3 = euc3;
		}
	}
	printf("there are %d inlier points\n", max_inlier_num);


	if (index_lst) { delete[]index_lst; index_lst = nullptr; }

};

