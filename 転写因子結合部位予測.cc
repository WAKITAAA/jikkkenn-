#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<cmath>
#include<random>
#include<algorithm>
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

//random配列のhit数を計算する関数
vector<vector<double>> HIT_random(vector<vector<double>> ozz){
     ifstream ist3("random_sequence");

    if(!ist3){
        cerr<<"Cannot open"<<endl;
        exit(1);
    }

    vector<string> promoter2;
    string s;
    vector <vector<double>> HIT(100); 
    int promoter_count=0;

    while(getline(ist3, s)){
        if(s[0]=='>'){
        promoter2.push_back(s);
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
    vector<vector<double>>ozz_MATalpha2=PFM(file_name[1]);
    vector<vector<double>>ozz_MCM1=PFM(file_name[2]);
    vector<vector<double>>ozz_MIG1=PFM(file_name[3]);
    vector<vector<double>>ozz_PHO4=PFM(file_name[4]);
    vector<vector<double>>ozz_RCS1=PFM(file_name[5]);
    vector<vector<double>>ozz_ROX1=PFM(file_name[6]);
    vector<vector<double>>ozz_TAF=PFM(file_name[7]);

    for(int i=0;i<4;i++){                        //ozz_MATa1の内容を出力してみる

        for(int j=0; j<ozz_MATa1[0].size(); j++){

            cout<<ozz_MATa1[i][j]<<" ";
        }
        cout<<endl;
    }
    

    vector<vector<double>>HIT_MATa1=HIT(ozz_MATa1);
    vector<vector<double>>HIT_MATalpha2=HIT(ozz_MATalpha2);
    vector<vector<double>>HIT_MCM1=HIT(ozz_MCM1);
    vector<vector<double>>HIT_MIG1=HIT(ozz_MIG1);
    vector<vector<double>>HIT_PHO4=HIT(ozz_PHO4);
    vector<vector<double>>HIT_RCS1=HIT(ozz_RCS1);
    vector<vector<double>>HIT_ROX1=HIT(ozz_ROX1);
    vector<vector<double>>HIT_TAF=HIT(ozz_TAF);

    cout<<HIT_MATa1[2][364]<<endl;
/*
    for(int i=0;i<8;i++){       //HIT_MATa1の内容を出力してみる

        for(int j=0; j<HIT_MATa1[0].size(); j++){

            cout<<HIT_MATa1[i][j]<<" ";
        }
        cout<<endl;
    }
    */

    cout<<"hello"<<endl;



    /*random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<> dist(0.0,1.0);

    vector<string> rannsuu_enki(100);

    for(int i = 0; i < 100; i++){      //１００本のランダム配列の生成
        rannsuu_enki[i]="";
        for(int j=0;j<500;j++){          //１本あたり、500塩基
            double ransuu=dist(mt);
            if((0.0<=ransuu)&&(ransuu<q_AT)){
                rannsuu_enki[i]+='A';
            }else if((q_AT<=ransuu)&&(ransuu<q_AT+q_CG)){
                rannsuu_enki[i]+='C';
            }else if((q_AT+q_CG<=ransuu)&&(ransuu<q_AT+2*q_CG)){
                rannsuu_enki[i]+='G';
            }else if((q_AT+2*q_CG<=ransuu)&&(ransuu<1.0)){
                rannsuu_enki[i]+='T';
            }
        }
    }

    ofstream fout("random_sequence");
    for(int i=0;i<100;i++){
        fout<<">random_sequence_"<<i+1<<endl;
        fout<<rannsuu_enki[i]<<endl;
    }
    fout.close();*/

    cout<<""<<endl;


     vector<vector<double>>HIT_random_MATa1=HIT_random(ozz_MATa1);//ランダム塩基配列１００本のhit数の計算　転写因子MAta1における
     vector<vector<double>>HIT_random_MATalpha2=HIT_random(ozz_MATalpha2);
     vector<vector<double>>HIT_random_MCM1=HIT_random(ozz_MCM1);
     vector<vector<double>>HIT_random_MIG1=HIT_random(ozz_MIG1);
     vector<vector<double>>HIT_random_PHO4=HIT_random(ozz_PHO4);
     vector<vector<double>>HIT_random_RCS1=HIT_random(ozz_RCS1);
     vector<vector<double>>HIT_random_ROX1=HIT_random(ozz_ROX1);
     vector<vector<double>>HIT_random_TAF=HIT_random(ozz_TAF);
/*
     for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_MATa1[0].size(); j++){

            cout<<HIT_random_MATa1[i][j]<<" ";
        }
        cout<<endl;
    }*/

    vector<double> hit_random_1dim_MATa1;      //ただ、hit数を詰め込んだ１次元配列
    vector<double> hit_random_1dim_MATalpha2;
    vector<double> hit_random_1dim_MCM1;
    vector<double> hit_random_1dim_MIG1;
    vector<double> hit_random_1dim_PHO4;
    vector<double> hit_random_1dim_RCS1;
    vector<double> hit_random_1dim_ROX1;
    vector<double> hit_random_1dim_TAF;

    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_MATa1[0].size(); j++){

            hit_random_1dim_MATa1.push_back(HIT_random_MATa1[i][j]);     //ランダム配列のhit数を１次元配列に変換
        }
    }
    sort(hit_random_1dim_MATa1.rbegin(), hit_random_1dim_MATa1.rend());   //ランダム配列のhit数を降順にソート


    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_MATalpha2[0].size(); j++){

            hit_random_1dim_MATalpha2.push_back(HIT_random_MATalpha2[i][j]);
        }
    }
    sort(hit_random_1dim_MATalpha2.rbegin(), hit_random_1dim_MATalpha2.rend());
    
    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_MCM1[0].size(); j++){

            hit_random_1dim_MCM1.push_back(HIT_random_MCM1[i][j]);
        }
    }
    sort(hit_random_1dim_MCM1.rbegin(), hit_random_1dim_MCM1.rend());

    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_MIG1[0].size(); j++){

            hit_random_1dim_MIG1.push_back(HIT_random_MIG1[i][j]);
        }
    }
    sort(hit_random_1dim_MIG1.rbegin(), hit_random_1dim_MIG1.rend());

    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_PHO4[0].size(); j++){

            hit_random_1dim_PHO4.push_back(HIT_random_PHO4[i][j]);
        }
    }
    sort(hit_random_1dim_PHO4.rbegin(), hit_random_1dim_PHO4.rend());

    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_RCS1[0].size(); j++){

            hit_random_1dim_RCS1.push_back(HIT_random_RCS1[i][j]);
        }
    }
    sort(hit_random_1dim_RCS1.rbegin(), hit_random_1dim_RCS1.rend());

    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_ROX1[0].size(); j++){

            hit_random_1dim_ROX1.push_back(HIT_random_ROX1[i][j]);
        }
    }
    sort(hit_random_1dim_ROX1.rbegin(), hit_random_1dim_ROX1.rend());

    for(int i=0;i<100;i++){

        for(int j=0; j<HIT_random_TAF[0].size(); j++){

            hit_random_1dim_TAF.push_back(HIT_random_TAF[i][j]);
        }
    }
    sort(hit_random_1dim_TAF.rbegin(), hit_random_1dim_TAF.rend());

    
    
    

    cout<<"転写因子MATa1のランダム配列のhit数の最大値:"<<hit_random_1dim_MATa1[0]<<endl;   //ランダム配列のhit数の最大値を出力 
    cout<<"ランダム配列の要素数:"<<hit_random_1dim_MATa1.size()<<endl;   //ランダム配列のhit数の要素数を出力
    
    //ｐ値を上位0.1パーとすると、上位４９個
    double p_value=0.00035;   //p値を0.2パーとする
    int index_MATa1=static_cast<int>(hit_random_1dim_MATa1.size()*p_value);   //p値に対応するインデックスを計算
    int index_MATalpha2=static_cast<int>(hit_random_1dim_MATalpha2.size()*p_value);
    int index_MCM1=static_cast<int>(hit_random_1dim_MCM1.size()*p_value);
    int index_MIG1=static_cast<int>(hit_random_1dim_MIG1.size()*p_value);
    int index_PHO4=static_cast<int>(hit_random_1dim_PHO4.size()*p_value);
    int index_RCS1=static_cast<int>(hit_random_1dim_RCS1.size()*p_value);
    int index_ROX1=static_cast<int>(hit_random_1dim_ROX1.size()*p_value);
    int index_TAF=static_cast<int>(hit_random_1dim_TAF.size()*p_value);     

    cout<<"転写因子MATa1のランダム配列のhit数の上位"<<p_value*100<<"パーの値:"<<hit_random_1dim_MATa1[index_MATa1]<<endl; 
    cout<<"転写因子MATalpha2のランダム配列のhit数の上位"<<p_value*100<<"パーの値:"<<hit_random_1dim_MATalpha2[index_MATalpha2]<<endl;
    cout<<"転写因子MCM1のランダム配列のhit数の上位"<<p_value*100<<"パーの値:"<<hit_random_1dim_MCM1[index_MCM1]<<endl;
    cout<<"転写因子MIG1のランダム配列のhit数の上位"<<p_value*100<<"パーの値:"<<hit_random_1dim_MIG1[index_MIG1]<<endl;
    cout<<"転写因子PHO4のランダム配列のhit数の上位"<<p_value*100<<"パーの値:"<<hit_random_1dim_PHO4[index_PHO4]<<endl;
    cout<<"転写因子RCS1のランダム配列のhit数の上位"<<p_value*100<<"パーの値:"<<hit_random_1dim_RCS1[index_RCS1]<<endl;
    cout<<"転写因子ROX1のランダム配列のhit数の上位"<<p_value*100<<"パーの   値:"<<hit_random_1dim_ROX1[index_ROX1]<<endl;
    cout<<"転写因子TAFのランダム配列のhit数の上位"<<p_value*100<<"パーの値:"<<hit_random_1dim_TAF[index_TAF]<<endl;     

    
    vector<string> promoter={"YDR270W","YJR048W","YDL227C","YBR093C","YPR065W","YFL026W","YKL209C","YLR377C"};   //プロモーター名を格納する配列
    cout<<endl;

    ifstream ist4("promoters");
    if(!ist4){
        cerr<<"Cannot open"<<endl;
        exit(1);
    }
   
    string s;
    vector<string> promoter_hairetu;

    while(getline(ist4, s)){
        if(s[0]!='>'){
            promoter_hairetu.push_back(s); 
        }
    }


    


    for(int i=0;i<8;i++){
        for(int j=0; j<HIT_MATa1[0].size(); j++){
            if(HIT_MATa1[i][j]>hit_random_1dim_MATa1[index_MATa1]){   //プロモーター配列のhit数がランダム配列の上位p%の値より大きい場合

                string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,10);
                cout<<"転写因子MATa1: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_MATa1[i][j]<<endl;   //プロモーター名、位置、hit数を出力

            }
        }
    }
    cout<<endl;
    
    for(int i=0;i<8;i++){
        for(int j=0; j<HIT_MATalpha2[0].size(); j++){
            if(HIT_MATalpha2[i][j]>hit_random_1dim_MATalpha2[index_MATalpha2]){
                 string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,10);
                cout<<"転写因子MATalpha2: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_MATalpha2[i][j]<<endl;
            }
        }
    }
    cout<<endl;

     for(int i=0;i<8;i++){
        for(int j=0; j<HIT_MCM1[0].size(); j++){
            if(HIT_MCM1[i][j]>hit_random_1dim_MCM1[index_MCM1]){
                 string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,16);
                cout<<"転写因子MCM1: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_MCM1[i][j]<<endl;
            }
        }
    }
    cout<<endl;

     for(int i=0;i<8;i++){
        for(int j=0; j<HIT_MIG1[0].size(); j++){
            if(HIT_MIG1[i][j]>hit_random_1dim_MIG1[index_MIG1]){
                string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,17);
                cout<<"転写因子MIG1: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_MIG1[i][j]<<endl;
            }
        }
    }
    cout<<endl;

     for(int i=0;i<8;i++){
        for(int j=0; j<HIT_PHO4[0].size(); j++){
            if(HIT_PHO4[i][j]>hit_random_1dim_PHO4[index_PHO4]){
                 string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,12);
                cout<<"転写因子PHO4: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_PHO4[i][j]<<endl;
            }
        }
    }
    cout<<endl;

     for(int i=0;i<8;i++){
        for(int j=0; j<HIT_RCS1[0].size(); j++){
            if(HIT_RCS1[i][j]>hit_random_1dim_RCS1[index_RCS1]){
                string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,13);
                cout<<"転写因子RCS1: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_RCS1[i][j]<<endl;
            }
        }
    }
    cout<<endl;

     for(int i=0;i<8;i++){
        for(int j=0; j<HIT_ROX1[0].size(); j++){
            if(HIT_ROX1[i][j]>hit_random_1dim_ROX1[index_ROX1]){
                string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,9);
                cout<<"転写因子ROX1: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_ROX1[i][j]<<endl;
            }
        }
    }
    cout<<endl;

    for(int i=0;i<8;i++){
        for(int j=0; j<HIT_TAF[0].size(); j++){
            if(HIT_TAF[i][j]>hit_random_1dim_TAF[index_TAF]){
                string hit_hairetu;
                hit_hairetu=promoter_hairetu[i].substr(j,15);
                cout<<"転写因子TAF: "<<promoter[i]<<" "<<j<<" "<<hit_hairetu<<" "<<HIT_TAF[i][j]<<endl;
            }
        }
    }




    return 0;
}
