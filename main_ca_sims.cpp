# include <cstdio>
# include <cstdlib>
# include <cmath>

#include <iostream> 
//#include <iterator>
//#include <algorithm>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <vector>
#include <experimental/filesystem>

using namespace std;

//** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
//** Victor Nicolai Friedhoff,  May 2023             **
//** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
//** Simulation programm for Ca2+ spikes, based on   **
//** a linear state chain w/ pos. feedback and       **
//** recovery from negative feedback.                **
//** See Falcke & Friedhoff (Chaos, 2018) and        **
//** Friedhoff et. al. (BPJ, 2023).                  **
//** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **

const string version_str="v1.00";

//***************** Classes ******************
class Parameters_sys {  //system parameters
    public:
        int max_state;      //if state max_state is reached, iteration ends.
        int runs;           //number of wanted succesfully finished iterations
        int runs_max;       //maximal number of iterations
        int max_time;       //if time > max_time, iteration gets aborted
        int save_intervall; //saves to file all save_intervall iterations
        double dt;          //time resolution/bin width
        double dt_base;     //back up of dt
};

class Parameters {          //model parameters
    public:
        double delta;       //down rate
        double sp;          //spatial parameter in pos. feedback
        double g;           //basis rate g(IP3)
        double Ka;          //for Hill positive feedback, Ka parameter in hill function

        //regular recovery
        double lambda;      //recovery rate from negative feedback
        double rr_n_rec;    //power of regular recovery term

        //recovery using hill function
        double rh_coeff;    //hill recovery, hill coefficient

        int Nt;             //total number of clusters
        int state;          //current state
        int n;              //exponent n in pos. feedback

        string pos_feedback_mode;
        string recovery_mode;
};


//**************** Functions ****************
//saves FPTs to file
void save_FPTs_vector(string out_path_file, vector<double> fpts, Parameters_sys sys_params, Parameters params) {
    ofstream write_out(out_path_file, ios::out | ios::app);

    //checking for opened stream
    if (!write_out.is_open()) {
        cerr << "Error opening file: " << out_path_file << endl;
        return;
    }

    //if file does not exist, add header built from file name
    string header, header_infos;
    header_infos = std::experimental::filesystem::path(out_path_file).stem();
    size_t dot_position = header_infos.find('.');
    header_infos = (dot_position != std::string::npos) ? header_infos.substr(dot_position + 1) : "";
    for (size_t found = header_infos.find('_'); found != std::string::npos; found = header_infos.find('_', found + 2)) {
        header_infos.replace(found, 1, ", ");
    }

    //tellp() returns the current position of the writing pointer in the file. If it's at position 0, the file is considered empty.
    if (write_out.tellp() == 0) {
        header = "# FPTs to states 0 up to max_state  // " + header_infos + "\n# ";
        for (int i = 0; i <= sys_params.max_state; i++) {
            header += to_string(i) + "     ";
        }
        header += "\n";
        write_out << header;
    }

    //dump FPTs into file
    int fpts_count = 0;
    for (double times : fpts) {
        write_out << fixed << setprecision(3) << times << " ";
        fpts_count += 1;
        if (fpts_count % (sys_params.max_state + 1) == 0) write_out << endl;
    }

    write_out.close();
}

//computes up rate based on positive feedback and recovery from negative feedback
double up_t_i(Parameters params, int state, double time){
    double recovery_term, out, pos_feedback_term;
    
    //********** recovery from negative feedback  **********
    //regular recovery from negative feedback, (1-e^(-lambda t))^n_rec
    recovery_term=pow((1-exp(-params.lambda*time)), params.rr_n_rec);
    
    //hill recovery from negative feedback, with x=(1-e^(-lambda t))^n_rec: (rh_c^n+1) * x^n/(rhc^n+x^n)
    if (params.recovery_mode=="hill"){
        recovery_term=(pow(params.rh_coeff,params.n)+1)*pow(recovery_term,params.n)/(pow(params.rh_coeff,params.n)+pow(recovery_term,params.n));
    }

    out = params.g * recovery_term * (params.Nt - state);

    //********** positive feedback  **********
    //regular positive feedback: (1+sp*state)^n
    pos_feedback_term = pow(1+params.sp*state,params.n);//

    //Hill positive feedback: 1 + sp * state^n / (Ka^n + state^n)
    if (params.pos_feedback_mode=="hill"){
        pos_feedback_term = 1+(params.sp*pow(state,params.n)/(pow(params.Ka,params.n)+pow(state,params.n)));
    }
    out = out * pos_feedback_term;

    return out;
}

//down rate
double down_t_i(double down_rate, int state){
    double out;

    out=down_rate*state;

    return out;
}

//for proper argument handling
bool is_num(const std::string& str) {
    std::istringstream iss(str);
    double value;
    iss >> std::noskipws >> value;
    return iss.eof() && !iss.fail();
}



//**************** Main Program ****************
int main(int argc, char* argv[])
{
    srand(time(NULL)); //chooses random seed on execution. Disable for making results reproduceable.

    //*********** variable handling ***********
    //**** parameter classes
    Parameters     params;
    Parameters_sys sys_params;

    //**** system variables and parameters
    bool   jump;
    double time, up_r, down_r, rate_check, p_up, p_do;     //time, rates, transition probabilities
    double time_out=100.0, prop_lim=0.15, succ_rate = 1;   //output time intervall, maximal transition probability before increasing time resoultion, success rate
    int    arg_num, check_folder, state, prev_state, highest_state=0, done_num=0;

    time_t time_now;
    std::chrono::time_point<std::chrono::system_clock> T_now;


    sys_params.runs_max = 15000;    //number of simulations (failed + successfull) before program terminates itself.
    sys_params.max_time = 4000;     //cancels simulation after max_time seconds if state Nt was not reached
    sys_params.save_intervall = 50; //saves to file after save_intervall simulations
    params.delta = 6.93;            //down rate

    //**** vectors for storing fpts
    vector<double> fpts = {};
    vector<double> fpts_buffer = {};

    //**** measuring time
    T_now    = std::chrono::system_clock::now();
    time_now = std::chrono::system_clock::to_time_t(T_now);
    auto time_start = std::chrono::high_resolution_clock::now();

    cout << "\n\033[1;33mMarkovian Calcium Spike Simulations, version " << version_str << ".\033[0m" << endl;
    cout <<   " Starting at " << std::ctime(&time_now) << endl;

    //*********** arguments ***********
    arg_num = 15;
    if (argc != arg_num-1) {
        cout << "\n\033[1;31m =======    Argument error     \033[1;31m======= \033[0m" << endl;
        cout << "   Need " << arg_num -1 << " arguments exactly, but found " << argc << "." << endl;
        cout << "   Positive Feedback should be 'exp' or 'hill', recovery mode 'normal' or 'hill'. Try" << endl;
        cout << "      ./caspikesims.out PosFdb   RecoveryMode Runs  dt  MaxState Lambda rr_nrec rh_coef Nt   g  sp  n Ka,   e.g.," << endl;
        cout << "      ./caspikesims.out exp      normal       20   1e-5   20     0.01    1.0      -1    40  3.0 1.0 3 -1" << endl;
        cout << "   rh_coef only applicable for Hill recovery, Ka only for Hill positive feedback." << endl;
        cout << "\n\033[1;31m =======    Ending program.    \033[1;31m======= \033[0m\n" << endl;
        return 0;
    }

    params.pos_feedback_mode   = argv[1];  //exp or hill, see up_t_i() function
    params.recovery_mode       = argv[2];  //normal or hill

    //check if other aguments are not strings.
    //attempt to convert to int and floats. if errors, its a string.
    for (int i = 3; i < argc; ++i) {
        if (!is_num(argv[i])) {
            cout << "   \033[1;31mError: \033[0m Some argument contains a string: " << argv[i] << ". Aborting." << endl;
            return 1;
        }
    }

    sys_params.runs      = atoi(argv[3]);  //number of simulations for given parameter set.
    sys_params.dt        = atof(argv[4]);  //initial time resolution
    sys_params.dt_base   = sys_params.dt;  //backup s.t. dt_base is used initially for each simulation iteration
    sys_params.max_state = atoi(argv[5]);  //simulation ends when state Nt or state max_state is reached
    params.lambda        = atof(argv[6]);  //recovery rate
    params.rr_n_rec      = atof(argv[7]);  //n_rec for regular recovery
    params.rh_coeff      = atof(argv[8]);  //Hill coefficient for hill recovery
    params.Nt            = atoi(argv[9]);  //total clusters in system
    params.g             = atof(argv[10]); //rate g(IP3)
    params.sp            = atof(argv[11]); //spatial parameter in positive feedback
    params.n             = atoi(argv[12]); //parameter n in positive feedback
    params.Ka            = atof(argv[13]); //only applicable for hill mode of positive feedback

    //optional. sets unused parameters to -1
    if ( params.recovery_mode      != "hill") params.rh_coeff = -1.0;
    if ( params.pos_feedback_mode  != "hill") params.Ka       = -1.0;
    
    cout << " called with \n    ";
    for(int i = 0; i < argc; ++i) cout << argv[i] << ' ';
    cout << endl << endl;

    cout << " Parameters: " << endl;
    cout << "    Pos Feedback & Recovery modes: " << params.pos_feedback_mode << ", " << params.recovery_mode << "\n" << endl;
    cout << "        runs = " << sys_params.runs     << ", \t     max_runs = " << sys_params.runs_max << endl;
    cout << "          dt = " << sys_params.dt       << ", \t    max_state = " << sys_params.max_state << endl;
    cout << "    max_time = " << sys_params.max_time << ", \t    down_rate = " << params.delta << endl;
    cout << "      Lambda = " << params.lambda << ", \t    rr_n_rec = " << params.rr_n_rec  <<  ", \t rh_coef = " << params.rh_coeff << endl;
    cout << "          Sp = " << params.sp     << ", \t           n = " << params.n         <<  ", \t      Nt = " << params.Nt << endl;
    cout << "          Ka = " << params.Ka     << ", \t           g = " << params.g         << endl << endl;


    //*********** checks and error handling  ***********
    if (params.Nt < sys_params.max_state) {
        cout << "   \033[1;31mError: \033[0mN_total < max_state. " << endl;
        cout << "          N_total > 10 is recommended. " << endl;    
        cout << "          => setting max_state to N_total." << endl;
        sys_params.max_state = params.Nt;
    }

    if (params.sp < 0 or params.lambda <= 0 or params.n < 0 or params.g <= 0 or params.Nt<=0 or params.rr_n_rec <0 or sys_params.runs <1 or sys_params.dt <=0) {
        cout << "   \033[1;31mError: \033[0mLambda, Sp, n, g, n_rec, Nt or time step dt is smaller than or equal to 0, or number of runs < 1. " << endl;
        return 1;
    }

    if (params.pos_feedback_mode == "hill" && params.Ka <= 0) {
        cout << "   \033[1;31mError: \033[0mKa is smaller than or equal to 0. " << endl;
        return 1;
    }

    if (params.recovery_mode == "hill" && params.rh_coeff <= 0) {
        cout << "   \033[1;31mError: \033[0mHill coefficient in recovery term is smaller than or equal to 0. " << endl;
        return 1;
    }

    if (params.recovery_mode != "normal" && params.recovery_mode != "hill") {
        cout << "   \033[1;31mError: \033[0mrecovery_mode \"" << params.recovery_mode << "\" not recognized. Use \"normal\" or \"hill\"." << endl;
        return 1;
    }

    if (params.pos_feedback_mode != "exp" && params.pos_feedback_mode != "hill") {
        cout << "   \033[1;31mError: \033[0mpos_feedback_mode \"" << params.pos_feedback_mode << "\" not recognized. Use \"exp\" or \"hill\"." << endl;
        return 1;
    }

    //******** string and folder handling ********
    string out_path, out_file, out_path_file, makefolder_str, recovery_str;
    string out_down_r, out_lambda, out_sp, out_g, out_dt, out_nrec, out_Ka;
    ostringstream ss_out_dt;

    out_nrec=to_string(params.rr_n_rec);
    out_nrec.erase( out_nrec.find_last_not_of('0') + 2, std::string::npos);

    out_path="results/PosFdbk_"+params.pos_feedback_mode+"__Recovery_"+params.recovery_mode+"__nrec=" +  out_nrec + "__Nt="+to_string(params.Nt);

    makefolder_str = "mkdir -p " + out_path;
    check_folder=system(makefolder_str.c_str());
    
    recovery_str="normal";
    if (params.recovery_mode=="hill") recovery_str="hill_rhcoeff="+to_string(params.rh_coeff);

    out_down_r  = to_string(params.delta);
    out_lambda  = to_string(params.lambda);
    out_sp      = to_string(params.sp);
    out_g       = to_string(params.g);
    out_Ka      = to_string(params.Ka);

    ss_out_dt.precision(3);
    ss_out_dt << scientific << sys_params.dt;
    out_dt = ss_out_dt.str();
   
    out_down_r.erase( out_down_r.find_last_not_of('0') + 2, std::string::npos);
    out_lambda.erase( out_lambda.find_last_not_of('0') + 1, std::string::npos);
    out_sp.erase(     out_sp.find_last_not_of('0') + 2,     std::string::npos);
    out_g.erase(      out_g.find_last_not_of('0') + 2,      std::string::npos);
    out_Ka.erase(     out_Ka.find_last_not_of('0') + 2,     std::string::npos);


    //******** handling file names ********
    if (params.pos_feedback_mode=="exp"){ //"_Up="+out_up_r+
        out_file="CaSim.PF="+params.pos_feedback_mode+"_Rec="+recovery_str+"_Down="+ out_down_r +"_n="+to_string(params.n) +"_nrec="+out_nrec+"_Lambda="+out_lambda+"_Nt="+to_string(params.Nt)  +"_Sp="+ out_sp +"_g="+out_g+"_MaxState="+to_string(sys_params.max_state)+"_dt="+out_dt+".txt";
        }
     else if (params.pos_feedback_mode=="hill"){
        out_file="CaSim.PF="+params.pos_feedback_mode+"_Rec="+recovery_str+"_Down="+out_down_r+"_n="+to_string(params.n)+"_nrec="+out_nrec+"_Lambda="+out_lambda+"_Nt="+to_string(params.Nt)+"_fmax="+out_sp+"_Ka="+ out_Ka +"_g="+out_g+"_MaxState="+to_string(sys_params.max_state)+"_dt="+out_dt+".txt";
    }
    out_path_file= out_path + "/" + out_file;

    cout << " Output filepath and filename:\n    " << out_path << endl;
    cout << "    " << out_file << "\n" <<endl;



    //***********
    //*********** actual simulation ***********
    //***********
    cout << "\033[1;34m ===== Starting simulations. ===== \033[0m" << endl;
    for (int run=0; run<sys_params.runs_max; run++) {
        if (run>0) succ_rate=(double)done_num/(double)run;

        if ((run>100) && succ_rate < 0.15){ //if after 100 iterations less than 15 have converged, abort program
            cout << "\n  Yielding few results. Only " << done_num << " of " << run << " runs finished." << endl;
            cout <<   "  That is a success rate of " << succ_rate * 100 << "\% < 15\% ->  Aborting simulations." <<endl;
            break;
        }

        //******* initialisation of simulation iteration
        sys_params.dt = sys_params.dt_base;
        time          = sys_params.dt;
        state         = 0;
        prev_state    = 0;
        highest_state = 0;
        fpts_buffer.push_back(0);

        cout << "   Run " << run+1 << "/"<< sys_params.runs << " (" << done_num << " done, "<<  round(succ_rate*100*100)/100 <<"\% success):     t \t i  \t  uprate  \t downrate \t dt" << endl;
        
        while (state < sys_params.max_state and time<sys_params.max_time) { //not converged and within timelimit
            up_r   = up_t_i(params, state, time);       //computes current up and down rates
            down_r = down_t_i(params.delta, state);
            
            p_up = up_r*sys_params.dt;                  //computes up and down transition/splitting probabilities
            if (p_up > prop_lim){ sys_params.dt=sys_params.dt/2; } //if too large, reduce bin width/time resolution dt
            
            p_up =   up_r*sys_params.dt;
            p_do = down_r*sys_params.dt;

            rate_check = up_r/(up_r+down_r);                //determines direction of jump
            jump = (p_up + p_do) > (float) rand()/RAND_MAX; //computes if jump occurs
            
            //cout << rate_check << ", " << p_up  << ", " << p_do<< ",  " <<  rand()/RAND_MAX << endl;
            
            //the jump
            if (jump==1){
                if ((float) rand()/RAND_MAX < rate_check) state+=1; 
                else if (state>0){ state-=1; } 
            } 
            
            if (state != prev_state){            //jumped?
                if (state > highest_state) {     //and new FPT?
                    fpts_buffer.push_back(time); //-> push FPT into buffer
                    //cout << "  reached state " << state << " at time " << time << endl;
                    highest_state=state;         //new highest state
                }
                if (state == sys_params.max_state) { done_num+=1; } //if finished, +1 succesfull iteration
            }
            
            //status echo
            if (fmod(time,time_out) < sys_params.dt or state == sys_params.max_state){
                cout << "                      \t       \t     " << time << " \t " << state << " \t " << up_r << " \t " << down_r << "  \t " << sys_params.dt << endl;
            }
            
            //handing over current values to next time bin
            prev_state=state;
            time+=sys_params.dt;
        }
        
        //has max_state been reached (finished iteration), i.e., is there an FPT for each state? add to FPTs & empty buffer
        if (fpts_buffer.size() == sys_params.max_state+1) { //+1 due to also counting state 0.
             fpts.insert(fpts.end(), fpts_buffer.begin(), fpts_buffer.end());
        }
        fpts_buffer={};

        //write current FPTs (sum of buffers) to file all save_intervall done iterations.
        //with save_intervall=1 and e.g. large g, Nt, Sp, ...  it might slow down MaxCluster due to many file operations
        if (done_num % sys_params.save_intervall == 0) {
            save_FPTs_vector(out_path_file, fpts, sys_params, params);
            fpts = {};
            cout << " ===== Saved. " << endl; 
        }

        cout << flush << endl;
        if (done_num == sys_params.runs) break; //end simulations if all iterations are done
    }

    save_FPTs_vector(out_path_file, fpts, sys_params, params); //save once more
    cout << " ===== Saved. " << endl;


    auto time_finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> time_elapsed = time_finish - time_start;
    cout << "\n   Successfully completed "  << done_num << " of " << sys_params.runs << " simulations." << endl;

    cout << setprecision(2);
    cout << "\n             Elapsed time: " << time_elapsed.count() << "s = "<< time_elapsed.count()/60 <<"m = " <<  time_elapsed.count()/3600 <<"h." <<  endl;
    cout <<   "   Average per Simulation: " << time_elapsed.count()/sys_params.runs << "s = "<< time_elapsed.count()/60/sys_params.runs<<"m.\n"<< endl;
    
    return 0;
}




