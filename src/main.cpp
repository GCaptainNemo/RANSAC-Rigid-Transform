#include "RANSAC_icp.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>



int main()
{
  std::vector<Eigen::Vector3d> src_pts(0);
  std::vector<Eigen::Vector3d> dst_pts(0);
  // inliers
  double src_lst[][3] = {
  	2.11224837, 7.19761707, 2.89912196, 6.62936141, 0.48177309,
       6.82718272, 0.83958316, 3.19141755, 6.59174479, 8.60479394,
       6.59092414, 0.22848416, 9.23657539, 9.51056172, 4.27600749,
       4.35313369, 4.07546213, 5.5735326 , 2.32060611, 9.25697726,
       7.39523616, 4.62254233, 8.53483148, 3.6943306 , 3.96160721,
       8.44185486, 7.53060809, 9.38447448, 3.74860614, 6.46442384,
       1.04557135, 9.5150682 , 1.45963384, 1.99798028, 3.19798554,
       1.29659021, 6.37363696, 2.42852721, 4.25923469, 7.5886749 ,
       2.86319074, 4.30121596, 6.06536736, 5.5502577 , 0.26595825,
       3.3686822 , 0.68834446, 5.160664  , 5.4664986 , 7.43293462,
       2.8448799 , 0.08065638, 5.59969693, 8.53764673, 1.97506394,
       7.9039781 , 7.10588289, 1.39005985, 6.07243229, 8.44892004
  };	
  Eigen::Matrix3d gt_R;
  gt_R << -0.55403972, -0.65322922, -0.51607323, 
           0.83244002, -0.44152023, -0.3348186 , 
           0.00914348,  0.61510282, -0.78839389;
  Eigen::Vector3d gt_T;
  gt_T << 8.61219902, 4.25520524, 1.90114084;
  Eigen::Isometry3d gt_iso = Eigen::Isometry3d::Identity();
  gt_iso.rotate(gt_R);
  gt_iso.pretranslate(gt_T);
  int inlier_N = sizeof(src_lst) / sizeof(double) / 3;
  for(int i = 0; i < inlier_N; ++i)
  {
      Eigen::Vector3d linshi;
      for(int j = 0; j < 3; ++j)
      {
         linshi(j, 0) = src_lst[i][j];
      }
      src_pts.push_back(linshi);
      Eigen::Vector3d dst_linshi = gt_iso * linshi;
      dst_pts.push_back(dst_linshi);
  }
  // add wrong match
  double src_noise[][3] = {2.01568926, 3.5204496 , 1.23081524, 1.86666879, 6.51238554,
       6.13785365, 5.87844325, 7.6050503 , 3.26221971, 7.84870335,
       8.34295926, 5.74325303, 2.19425418, 0.81883556, 6.02722213,
       5.39859708, 9.19333587, 4.17618172, 3.85877432, 1.45460673,
       3.943413  , 1.01936432, 2.02737386, 6.7234084 , 9.47051896,
       8.41808683, 6.28187952, 8.35167143, 2.0021455 , 0.98637587};
  double dst_noise[][3] = {
       8.90549342, 6.20835234, 7.78359823, 2.33693187, 9.71638037,
       2.59899435, 5.5954664 , 4.39870832, 5.70814874, 3.92676507,
       7.15993713, 1.94932906, 6.40784067, 3.18898044, 8.71665891,
       9.86751143, 1.9118042 , 0.90111544, 1.1766824 , 8.49101835,
       7.37003675, 6.96288139, 5.20249482, 0.78760861, 6.72701428,
       0.67754361, 3.19502269, 2.54959072, 6.33217001, 7.37481951};
  int outlier_N = sizeof(src_noise) / sizeof(double) / 3;
  for(int i = 0; i < outlier_N; ++i)
  {
      Eigen::Vector3d src_linshi;
      Eigen::Vector3d dst_linshi;

      for(int j = 0; j < 3; ++j)
      {
         src_linshi(j, 0) = src_noise[i][j];
         dst_linshi(j, 0) = dst_noise[i][j];
      }
      src_pts.push_back(src_linshi);
      dst_pts.push_back(dst_linshi);
  }
  RansacRigidTransform ransac_obj;
  std::vector<uchar> status(0);
  ransac_obj.start_RANSAC(src_pts, dst_pts, status);
  std::cout << "final result euc3 = " << std::endl;
  std::cout << ransac_obj.res_euc3.matrix() << std::endl;
  std::cout << "delta euc3 = " << std::endl;
  std::cout << ransac_obj.res_euc3.matrix() -  gt_iso.matrix()<< std::endl;

  return 0;
}
