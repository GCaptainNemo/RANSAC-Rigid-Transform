#include "RANSAC_icp.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>



int main()
{
  std::vector<Eigen::Vector3d> src_pts(0);
  std::vector<Eigen::Vector3d> dst_pts(0);
  double src_lst[10][3] = {
  	2.11224837, 7.19761707, 2.89912196, 6.62936141, 0.48177309,
       6.82718272, 0.83958316, 3.19141755, 6.59174479, 8.60479394,
       6.59092414, 0.22848416, 9.23657539, 9.51056172, 4.27600749,
       4.35313369, 4.07546213, 5.5735326 , 2.32060611, 9.25697726,
       7.39523616, 4.62254233, 8.53483148, 3.6943306 , 3.96160721,
       8.44185486, 7.53060809, 9.38447448, 3.74860614, 6.46442384
  };	
  Eigen::Matrix3d gt_R;
  gt_R << -0.55403972, -0.65322922, -0.51607323, 
           0.83244002, -0.44152023, -0.3348186 , 
           0.00914348,  0.61510282, -0.78839389;
  Eigen::Vector3d gt_T;
  gt_T << 8.61219902, 4.25520524, 1.90114084;
  for(int i=0; i < 10; ++i)
  {
      Eigen::Vector3d linshi;
      for(int j =0; j < 3; ++j)
      {
         linshi(j, 0) = src_lst[i][j];
      }
      src_pts.push_back(linshi);
      Eigen::Vector3d dst_linshi = gt_R * linshi + gt_T;
      dst_pts.push_back(dst_linshi);
  }
  RansacRigidTransform ransac_obj;
  std::vector<uchar> status(0);
  ransac_obj.start_RANSAC(src_pts, dst_pts, status);
  std::cout << "final result euc3 = " << std::endl;
  std::cout << ransac_obj.res_euc3.matrix() << std::endl;

  return 0;
}
