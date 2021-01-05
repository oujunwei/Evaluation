#ifndef CONFIG_H
#define CONFIG_H

#define PI 3.141592653589793238462643383279502884197169399375105

using namespace std;

vector<double> probreference;

int nt=10;//环境变化强度
int run=5;    //独立运行次数
int predict=100; //预测代数
int metrics=5; //指标数
int testNumber=14;//测试函数个数
int precise = 3;  //精确度
//char* instances[]  = {"FDA1","FDA2","FDA3","FDA4","DMOP1","DMOP2","DMOP3","DMOPA","DMOPB","DMOPC","DMOPD","DMOPE","JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9","FDA4","DMOPF"}; // names of test instances
	
//char* instances[]  = {"FDA1","FDA2","FDA3","FDA4","FDA5","DMOP1","DMOP2","DMOP3","JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9"}; // names of test instances
char* instances[]  = {"DF1","DF2","DF3","DF4","DF5","DF6","DF7","DF8","DF9","DF10","DF11","DF12","DF13","DF14"};

#endif // !CONFIG_H