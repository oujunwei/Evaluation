#include "evaluate.h"
void average();

void stateAverage();

void readata(char *filename, vector< vector<double>> &pf);

using namespace std;

int main()
{
	average();
	return 0;
}


void readata(char *filename, vector< vector<double>> &pf){
	std::fstream fin;
	int line=0;
	char str[100]=" ";
	fin.open(filename,std::ios::in);
	pf.clear();
	if(fin.is_open())
	{
		const char* pch2;		
		std::string str;
		while(!fin.eof())
		{
			vector<double> data;
			std::getline(fin,str,'\n');
			pch2 = str.c_str();
			const char *temp=" ";
			char *p;
			char* ccha = nullptr; 
			ccha= const_cast<char*>(pch2); 
			
			p=strtok(ccha,temp);
			while (p)
			{				
				data.push_back(atof(p));
				p=strtok(NULL,temp);
			}
			line++;
			pf.push_back(data);
		}
	} //end if

	else
		std::cout<<"failed to open "<<filename<<endl;
	fin.close();
}


void average(){

	int ccc=0;
	std::ofstream pf;
	std::ofstream sta;

	while(ccc<run)
	{
		evaluate eva;
		
		vector<double> gd;
		char filename1[1024];
		char filename2[1024];
		char    strTestInstance[256];
		for (int i=0;i<testNumber;i++)
		{			
			if(instances[i]=="JY1"||instances[i]=="JY2"||instances[i]=="JY3"||instances[i]=="JY4"||instances[i]=="JY5"
				||instances[i]=="JY6"||instances[i]=="JY7"||instances[i]=="JY8"||instances[i]=="JY9"||instances[i]=="DMOP3"||
				instances[i]=="FDA1"||instances[i]=="FDA2"||instances[i]=="FDA3"||instances[i]=="DMOP1"||instances[i]=="DMOP2"||
				instances[i]=="DMOPA"||instances[i]=="DMOPB"||instances[i]=="DMOPC"||instances[i]=="DMOPD"||instances[i]=="DMOPE"||
				instances[i]=="DF1" || instances[i]=="DF2"||instances[i]=="DF3" || instances[i]=="DF4" || instances[i]=="DF5"||
				instances[i]=="DF6" || instances[i]=="DF7"||instances[i]=="DF8" || instances[i]=="DF9" ){
				char pfname[1024];
				sprintf(pfname,"evaluate/data/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);			
				
				for(int gen=1; gen<=predict;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i]);
					//sprintf(filename1,"PF/pf_%s.dat",strTestInstance);
					//eva.loadpfront(filename1,eva.pf,2);	//真实的PF放到文件中	
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点

					sprintf(filename2,"PF/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);
					eva.loadpfront(filename2,eva.p,2);

					//pf<<setprecision(7)<<setiosflags(ios::scientific);
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					pf<<eva.indicator_hvd(eva.p,gen,instances[i],nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}				
			else if(instances[i]=="DMOPF" ||instances[i]=="FDA4" ||instances[i]=="FDA5"|| instances[i]=="DF10"
				 || instances[i]=="DF11"||instances[i]=="DF12" || instances[i]=="DF13" ||instances[i]=="DF14")
			{
				char pfname[1024];
				double ave=0.0;
				int gen=0;
				sprintf(pfname,"evaluate/data/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);

				for(gen=1; gen<=predict;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i]);
					
					
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点

					sprintf(filename2,"PF/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);										
					eva.loadpfront(filename2,eva.p,3);

					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					pf<<eva.indicator_hvd(eva.p,gen,instances[i],nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;			
					eva.p.clear();
					eva.pf.clear();
				}			
				pf.close();
			}
			
		}
		ccc++;
	}
	

	//求样本平均值和样本方差	
	evaluate eva;
	for (int kkk = 0 ; kkk < testNumber; kkk++)
	{
		char pfname[1024];
		
		double ave=0.0;
		int t=0;
		sprintf(pfname,"evaluate/avg/%s.dat",instances[kkk]);
		
		pf.open(pfname,ios::trunc);

		char    strTestInstance[256];
		char    strTestInstanceAv[256];

		cout<<instances[kkk]<<endl;

		pf<<setprecision(5)<<setiosflags(ios::scientific);
		for(ccc=0;ccc<run;ccc++){     //每次独立运行的求平均
			sprintf(strTestInstance,"evaluate/data/%s_%d.dat",instances[kkk],ccc);
			eva.loadpfront(strTestInstance,eva.pf,metrics);
			eva.indicator_AVG(eva.pf,eva.in_avg);
			for(unsigned  int iii=0;iii<eva.in_avg.size();iii++){
				pf<<eva.in_avg[iii]<<"      ";
			}
			pf<<"\n";	
		}
		pf.close();

		//独立run次的平均和方差
		pf.open(pfname,ios::app);
		pf<<setprecision(5)<<setiosflags(ios::scientific);
		sprintf(strTestInstanceAv,"evaluate/avg/%s.dat",instances[kkk]);
		eva.loadpfront(strTestInstanceAv,eva.pf,metrics);
		eva.indicator_AVG(eva.pf,eva.in_avg);
		for(unsigned  int iii=0;iii<eva.in_avg.size();iii++){
				pf<<eva.in_avg[iii]<<"      ";
		}
		pf<<"\n";
		unsigned  int size1=eva.pf.size();
		unsigned  int size2=eva.pf[0].size();
		vector<double> values;
		for(unsigned  int iii=0;iii<size2;iii++){
			for(unsigned ccc=0;ccc<size1;ccc++){
				values.push_back(eva.pf[ccc][iii]);
			}
			pf<<eva.STD(values)<<"      ";
			values.clear();
		}
		pf<<"\n";
		pf.close();

		//求平均IGD
		eva.GetAvgIGD(instances[kkk],run,predict);
		
	}
	eva.Statistics(instances,testNumber,metrics,precise,run);//统计数据
}

void stateAverage(){
	int ccc=0;
	int nt=10;          //环境变化强度
	int run=5;          //独立运行次数
	int predict=120;    //预测代数
	int metrics=4;      //指标数
	int testNumber=2;   //测试函数个数

	std::ofstream pf;
	std::ofstream sta;

	char* instances[]  = {"FDA1","DMOP1","FDA2","FDA3","FDA4","DMOP2","DMOP3","DMOPA","DMOPB","DMOPC","DMOPD","DMOPE","DMOPF","JY1","JY2","JY3","JY4","JY5"}; // names of test instances
	
	//char* instances[]  = {"JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9"}; // names of test instances
	
	while(ccc<run)
	{
		evaluate eva;
			
		vector<double> gd;
		
		char filename2[1024];
		char    strTestInstance[256];
		for (int i=0; i< testNumber;i++)
		{
			//二维
			if (instances[i]=="JY1"||instances[i]=="JY2"||instances[i]=="JY3"||instances[i]=="JY4"||instances[i]=="JY5"
				||instances[i]=="JY6"||instances[i]=="JY7"||instances[i]=="JY8"||instances[i]=="JY9"||
				instances[i]=="FDA1" ||instances[i]=="FDA2"||instances[i]=="FDA3" ||instances[i]=="DMOP1"||instances[i]=="DMOP2"
				||instances[i]=="DMOP3"||instances[i]=="DMOPA" ||instances[i]=="DMOPB"||instances[i]=="DMOPC"
				||instances[i]=="DMOPD"||instances[i]=="DMOPE" )
			{
				char pfname[1024];
				sprintf(pfname,"state/1st/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);

				double ave=0.0,fangcha = 0.0;
				int gen=0;
				for(gen=1; gen<=20;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i]);					
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点
					
					sprintf(filename2,"P/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);
					eva.loadpfront(filename2,eva.p,2);
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               "<<"\n";
					//<<eva.indicator_hvd(eva.p,t,instances[i],nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;

					eva.p.clear();
					eva.pf.clear();
				}				
				pf.close();
			}
			if (instances[i]=="JY1"||instances[i]=="JY2"||instances[i]=="JY3"||instances[i]=="JY4"||instances[i]=="JY5"
				||instances[i]=="JY6"||instances[i]=="JY7"||instances[i]=="JY8"||instances[i]=="JY9"||
				instances[i]=="FDA1" ||instances[i]=="FDA2"||instances[i]=="FDA3" ||instances[i]=="DMOP1"||instances[i]=="DMOP2"
				||instances[i]=="DMOP3"||instances[i]=="DMOPA" ||instances[i]=="DMOPB"||instances[i]=="DMOPC"
				||instances[i]=="DMOPD"||instances[i]=="DMOPE" )
			{
				char pfname[1024];
				sprintf(pfname,"state/2st/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);

				double ave=0.0,fangcha = 0.0;
				int gen=0;
				for(gen=21; gen<=40;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i]);					
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点
					
					sprintf(filename2,"P/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);
					eva.loadpfront(filename2,eva.p,2);
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               "<<"\n";
					//<<eva.indicator_hvd(eva.p,t,instances[i],nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;

					eva.p.clear();
					eva.pf.clear();
				}				
				pf.close();
			}
			if (instances[i]=="JY1"||instances[i]=="JY2"||instances[i]=="JY3"||instances[i]=="JY4"||instances[i]=="JY5"
				||instances[i]=="JY6"||instances[i]=="JY7"||instances[i]=="JY8"||instances[i]=="JY9"||
				instances[i]=="FDA1" ||instances[i]=="FDA2"||instances[i]=="FDA3" ||instances[i]=="DMOP1"||instances[i]=="DMOP2"
				||instances[i]=="DMOP3"||instances[i]=="DMOPA" ||instances[i]=="DMOPB"||instances[i]=="DMOPC"
				||instances[i]=="DMOPD"||instances[i]=="DMOPE" )
			{
				char pfname[1024];
				sprintf(pfname,"state/3st/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);

				double ave=0.0,fangcha = 0.0;
				int gen=0;
				for(gen=41; gen<=80;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i]);					
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点
					
					sprintf(filename2,"P/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);
					eva.loadpfront(filename2,eva.p,2);
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               "<<"\n";
					//<<eva.indicator_hvd(eva.p,t,instances[i],nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;

					eva.p.clear();
					eva.pf.clear();
				}				
				pf.close();
			}
			//三维
			if (instances[i]=="FDA4"||instances[i]=="DMOPF")
			{
				char pfname[1024];
				sprintf(pfname,"state/1st/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);

				double ave=0.0,fangcha = 0.0;
				int t=0;
				for(t=1; t<=40;t++)
				{
					sprintf(strTestInstance,"%s",instances[i]);				
					eva.getPOF(instances[i],nt,t,eva.pf);  //取真实的点

					sprintf(filename2,"P/pf_%s_%d_%d.dat",strTestInstance,ccc,t);										
					eva.loadpfront(filename2,eva.p,3);

					//pf<<setprecision(7)<<setiosflags(ios::scientific);
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					//pf<<eva.indicator_SP(eva.p)<<"               ";
					//pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					//<<eva.indicator_hvd(eva.p,t,instances[i],nt)<<"\n";
					pf<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  t = "<<t<<"  "<<endl;
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}
			if (instances[i]=="FDA4"||instances[i]=="DMOPF")
			{
				char pfname[1024];
				sprintf(pfname,"state/2st/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);

				double ave=0.0,fangcha = 0.0;
				int t=0;
				for(t=41; t<=80;t++)
				{
					sprintf(strTestInstance,"%s",instances[i]);				
					eva.getPOF(instances[i],nt,t,eva.pf);  //取真实的点

					sprintf(filename2,"P/pf_%s_%d_%d.dat",strTestInstance,ccc,t);										
					eva.loadpfront(filename2,eva.p,3);

					//pf<<setprecision(7)<<setiosflags(ios::scientific);
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					//pf<<eva.indicator_SP(eva.p)<<"               ";
					//pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					//<<eva.indicator_hvd(eva.p,t,instances[i],nt)<<"\n";
					pf<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  t = "<<t<<"  "<<endl;
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}		
			if (instances[i]=="FDA4"||instances[i]=="DMOPF")
			{
				char pfname[1024];
				sprintf(pfname,"state/3st/%s_%d.dat",instances[i],ccc);

				pf.open(pfname,ios::trunc);

				double ave=0.0,fangcha = 0.0;
				int t=0;
				for(t=81; t<=120;t++)
				{
					sprintf(strTestInstance,"%s",instances[i]);				
					eva.getPOF(instances[i],nt,t,eva.pf);  //取真实的点

					sprintf(filename2,"P/pf_%s_%d_%d.dat",strTestInstance,ccc,t);										
					eva.loadpfront(filename2,eva.p,3);

					//pf<<setprecision(7)<<setiosflags(ios::scientific);
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					//pf<<eva.indicator_SP(eva.p)<<"               ";
					//pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					//<<eva.indicator_hvd(eva.p,t,instances[i],nt)<<"\n";
					pf<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  t = "<<t<<"  "<<endl;
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}		
		}
		ccc++;
	}
	
	//1st
	evaluate eva1;
	for (int kkk = 0 ; kkk < testNumber; kkk++)
	{
		char pfname[1024];
		
		double ave=0.0;
		int t=0;
		sprintf(pfname,"state/1st/avg/%s.dat",instances[kkk]);
		
		pf.open(pfname,ios::trunc);

		char    strTestInstance[256];
		char    strTestInstanceAv[256];

		cout<<instances[kkk]<<endl;

		pf<<setprecision(3)<<setiosflags(ios::scientific);
		for(ccc=0;ccc<run;ccc++){     //每次独立运行的求平均
			sprintf(strTestInstance,"state/1st/%s_%d.dat",instances[kkk],ccc);
			eva1.loadpfront(strTestInstance,eva1.pf,metrics);
			eva1.indicator_AVG(eva1.pf,eva1.in_avg);
			for(unsigned int iii=0;iii<eva1.in_avg.size();iii++){
				pf<<eva1.in_avg[iii]<<"      ";
			}
			pf<<"\n";	
		}
		pf.close();

		//独立run次的平均和方差
		pf.open(pfname,ios::app);
		pf<<setprecision(3)<<setiosflags(ios::scientific);
		sprintf(strTestInstanceAv,"state/1st/avg/%s.dat",instances[kkk]);
		eva1.loadpfront(strTestInstanceAv,eva1.pf,metrics);
		eva1.indicator_AVG(eva1.pf,eva1.in_avg);
		for(unsigned int iii=0;iii<eva1.in_avg.size();iii++){
				pf<<eva1.in_avg[iii]<<"      ";
		}
		pf<<"\n";
		int size1=eva1.pf.size();
		int size2=eva1.pf[0].size();
		vector<double> values;
		for(int iii=0;iii<size2;iii++){
			for(ccc=0;ccc<size1;ccc++){
				values.push_back(eva1.pf[ccc][iii]);
			}
			pf<<eva1.STD(values)<<"      ";
			values.clear();
		}
		pf<<"\n";
		pf.close();
	}

	//2st
	evaluate eva2;
	for (int kkk = 0 ; kkk < testNumber; kkk++)
	{
		char pfname[1024];
		
		double ave=0.0;
		int t=0;
		sprintf(pfname,"state/2st/avg/%s.dat",instances[kkk]);
		
		pf.open(pfname,ios::trunc);

		char    strTestInstance[256];
		char    strTestInstanceAv[256];

		cout<<instances[kkk]<<endl;

		pf<<setprecision(3)<<setiosflags(ios::scientific);
		for(ccc=0;ccc<run;ccc++){     //每次独立运行的求平均
			sprintf(strTestInstance,"state/2st/%s_%d.dat",instances[kkk],ccc);
			eva2.loadpfront(strTestInstance,eva2.pf,metrics);
			eva2.indicator_AVG(eva2.pf,eva2.in_avg);
			for(unsigned int iii=0;iii<eva2.in_avg.size();iii++){
				pf<<eva2.in_avg[iii]<<"      ";
			}
			pf<<"\n";	
		}
		pf.close();

		//独立run次的平均和方差
		pf.open(pfname,ios::app);
		pf<<setprecision(3)<<setiosflags(ios::scientific);
		sprintf(strTestInstanceAv,"state/2st/avg/%s.dat",instances[kkk]);
		eva2.loadpfront(strTestInstanceAv,eva2.pf,metrics);
		eva2.indicator_AVG(eva2.pf,eva2.in_avg);
		for(unsigned int iii=0;iii<eva2.in_avg.size();iii++){
				pf<<eva2.in_avg[iii]<<"      ";
		}
		pf<<"\n";
		int size1=eva2.pf.size();
		int size2=eva2.pf[0].size();
		vector<double> values;
		for(int iii=0;iii<size2;iii++){
			for(ccc=0;ccc<size1;ccc++){
				values.push_back(eva2.pf[ccc][iii]);
			}
			pf<<eva2.STD(values)<<"      ";
			values.clear();
		}
		pf<<"\n";
		pf.close();
	}

	//3st
	evaluate eva3;
	for (int kkk = 0 ; kkk < testNumber; kkk++)
	{
		char pfname[1024];
		
		double ave=0.0;
		int t=0;
		sprintf(pfname,"state/3st/avg/%s.dat",instances[kkk]);
		
		pf.open(pfname,ios::trunc);

		char    strTestInstance[256];
		char    strTestInstanceAv[256];

		cout<<instances[kkk]<<endl;

		pf<<setprecision(3)<<setiosflags(ios::scientific);
		for(ccc=0;ccc<run;ccc++){     //每次独立运行的求平均
			sprintf(strTestInstance,"state/3st/%s_%d.dat",instances[kkk],ccc);
			eva3.loadpfront(strTestInstance,eva3.pf,metrics);
			eva3.indicator_AVG(eva3.pf,eva3.in_avg);
			for(unsigned int iii=0;iii<eva3.in_avg.size();iii++){
				pf<<eva3.in_avg[iii]<<"      ";
			}
			pf<<"\n";	
		}
		pf.close();

		//独立run次的平均和方差
		pf.open(pfname,ios::app);
		pf<<setprecision(3)<<setiosflags(ios::scientific);
		sprintf(strTestInstanceAv,"state/3st/avg/%s.dat",instances[kkk]);
		eva3.loadpfront(strTestInstanceAv,eva3.pf,metrics);
		eva3.indicator_AVG(eva3.pf,eva3.in_avg);
		for(unsigned  int iii=0;iii<eva3.in_avg.size();iii++){
				pf<<eva3.in_avg[iii]<<"      ";
		}
		pf<<"\n";
		int size1=eva3.pf.size();
		int size2=eva3.pf[0].size();
		vector<double> values;
		for(int iii=0;iii<size2;iii++){
			for(ccc=0;ccc<size1;ccc++){
				values.push_back(eva3.pf[ccc][iii]);
			}
			pf<<eva3.STD(values)<<"      ";
			values.clear();
		}
		pf<<"\n";
		pf.close();
	}
}

