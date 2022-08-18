# SPAC
Code for "SPAC: Sparse Sensor Placement Based Adaptive Control for High Precision Fuselage Assembly" 
CODE DESCRIPTION 

*The main function is "AdaptiveAndSparseSensing.m". Please use Matlab (we use R2021a, any recent version should work) to run it.

*We use dummy data to replace the Boeing data. Therefore, running this code does not reproduce the results in the paper because the data is not the original real Boeing data. If you wish to use your own data, there are four types of data to be prepared:
1) stiffness matrix of the designed shape(K0)
2) stiffness matrix of the true part(K_t)
3) initial deformation vector (d0)
4) gravity force vector (g)
All of them can be directly exported from FEA software. A detailed description can be found in the code.


*Licensing information (BSD License)
Link to code (The code will be deposited to the authors' webpage if the manuscript is accepted. It will also be provided as supplementary materials accompanying this paper)


* Note: to run the code, cvx software is needed (http://cvxr.com/). 


REPRODUCIBILITY
Running this code does not reproduce the results in the paper because the data is not the original real Boeing data.


CONTACT INFORMATION 

Jianjun (Jan) Shi
The Carolyn J. Stewart Chair and Professor
H. Milton Stewart School of Industrial and Systems Engineering
765 Ferst Drive, Groseclose Building, Room 109
Georgia Institute of Technology
Atlanta, Georgia 30332-0205
jianjun.shi@isye.gatech.edu 
Phone:  404-385-3488, Fax:  404 894 2301
http://www2.isye.gatech.edu/~jshi33

Xiaowei Yue
Assistant Professor
Grado Department of Industrial & Systems Engineering
Virginia Tech, 113 Durham Hall (MC 0118)
Blacksburg, VA 24061
E-mail: xwy@vt.edu; Phone: 540-231-9081
Website: http://ise.vt.edu/yue

Shancong Mou
H. Milton Stewart School of Industrial and Systems Engineering
765 Ferst Drive, Groseclose Building, Room 109
Georgia Institute of Technology
Atlanta, Georgia 30332-0205
E-mail: shancong.mou@gatech.edu; Phone: 404-200-5070
