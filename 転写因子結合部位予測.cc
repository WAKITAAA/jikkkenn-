#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<cmath>
using namespace std;

const double total_bg=24314210.0;
const double bg_at=7519429.0;
const double bg_cg=4637676.0;
const double q_AT=bg_at/total_bg;
const double q_CG=bg_cg/total_bg;


vector <string> file_name={"MATa1","MATalpha2","MCM1","MIG1","PHO4","RCS1","ROX1","TAF"};

   

vector<vector<double>> PFM(string File_name){

    ifstream ist(File_name);
    if(!ist){
        cerr<<"Cannot open"<<endl;
        exit(1);
    }
            
    string s;
    vector <vector<double>> atgc_count(4,vector<double>(20,1.0)); //１で初期化しとく

    int length;
    while(getline(ist, s)){
        length=s.size();
        for(int i=0;i<s.size();i++){
            if(s[i]=='A'){
            atgc_count[0][i]++;
            }else if(s[i]=='C'){
            atgc_count[1][i]++;
            }else if(s[i]=='G'){
                atgc_count[2][i]++;
            }else if(s[i]=='T'){
                atgc_count[3][i]++;
            }
        }
    }

      for(int i=0;i<4;i++){

        for(int j=0; j<atgc_count[0].size(); j++){

            cout<<atgc_count[i][j]<<" ";
        }
        cout<<endl;
    }

 
    vector <vector<double>> ozz(4,vector<double>(length,0));
    
    for(int j=0;j<length;j++){
        int total_bp=0;

        for(int i=0;i<4;i++){
            total_bp=total_bp+atgc_count[i][j];
        }

        for(int i=0;i<4;i++){

        if((i==0)||(i==3)){
       

        ozz[i][j]=log((atgc_count[i][j]/ total_bp)/q_AT);
        }else{
        ozz[i][j]=log((atgc_count[i][j]/ total_bp)/q_CG);
        }
    }
        }

    return ozz;
 }


vector<vector<double>> HIT(vector<vector<double>> ozz){
     ifstream ist2("promoters");

    if(!ist2){
        cerr<<"Cannot open"<<endl;
        exit(1);
    }

    vector<string> promoter;
    string s;
    vector <vector<double>> HIT(8); 
    int promoter_count=0;

    while(getline(ist2, s)){
        if(s[0]=='>'){
        promoter.push_back(s);
        promoter_count++;

    }else if(s[0]!='>'){

        for(int i=0;i<s.size()-ozz[0].size()+1;i++){
            double hit=0;
            for(int k=i;k<i+ozz[0].size();k++){
                if(s[k]=='A'){
                    hit=hit+ozz[0][k-i];
                }else if(s[k]=='C'){
                    hit=hit+ozz[1][k-i];
                }else if(s[k]=='G'){
                    hit=hit+ozz[2][k-i];
                }else if(s[k]=='T'){
                    hit=hit+ozz[3][k-i];
                }

            }
            HIT[promoter_count-1].push_back(hit);
        }
    }
    }
    return HIT;
}





int main(void){

    vector<vector<double>>ozz_MATa1=PFM(file_name[0]);

    for(int i=0;i<4;i++){

        for(int j=0; j<ozz_MATa1[0].size(); j++){

            cout<<ozz_MATa1[i][j]<<" ";
        }
        cout<<endl;
    }

    vector<vector<double>>HIT_MATa1=HIT(ozz_MATa1);

    cout<<HIT_MATa1[2][364]<<endl;

    for(int i=0;i<8;i++){

        for(int j=0; j<HIT_MATa1[0].size(); j++){

            cout<<HIT_MATa1[i][j]<<" ";
        }
        cout<<endl;
    }




    return 0;
}
