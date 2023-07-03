#include "eigen3/Eigen/Dense"
#include "RicVolumeSet.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <omp.h>
#include <list>
using namespace std;
const double K_PI = 3.1415926535897;

static list<complex<double>> pad_list(list<complex<double>> x){
	const int N = x.size();
	int M = 1;
	while(M < N){
		M *= 2;
	}
	const int diff = M - N;
	if(diff){
		for(int i = 0; i < diff; ++i){
			complex<double> zero = 0.0;
			x.push_back(zero);
		}
	}
	return x;
}

static list<complex<double>> radix2_FFT(list<complex<double>> x){
	
	int N = x.size();

	if(N > 1){

		list<complex<double>> x_even(N/2);
		list<complex<double>> x_odd(N/2);
		list<complex<double>>::iterator x_iter = x.begin();
		list<complex<double>>::iterator x_even_iter = x_even.begin();
		list<complex<double>>::iterator x_odd_iter = x_odd.begin();

		while(x_iter != x.end()){
			*x_even_iter++ = *x_iter++;
			*x_odd_iter++ = *x_iter++;
		}
		
		x_even = radix2_FFT(x_even);
		x_odd = radix2_FFT(x_odd);
		complex<double> omega = polar(1.0, -2.0*K_PI/N);
		
		x_odd_iter = x_odd.begin();
		int k = 0;
		list<complex<double>>::iterator x_iter_stride = x.begin();
		while(x_odd_iter != x_odd.end()){
			*x_odd_iter++ *= pow(omega, k++);
			++x_iter_stride;
		}		

		x_odd_iter = x_odd.begin();
		x_even_iter = x_even.begin();
		x_iter = x.begin();

		while(x_even_iter != x_even.end()){
			const complex<double> x_even_value = *x_even_iter++;
			const complex<double> x_odd_value = *x_odd_iter++;
			*x_iter++ = x_even_value + x_odd_value;
			*x_iter_stride++ = x_even_value - x_odd_value;
		}
	}
	return x;
}

static void print_help(){
	 cout << "reho_csf -i <input volume name> -mask <mask volume name> -o <output file base name>\n \
	  -t <time per image>\n";
	 cout << " optional arguments: -bins <number of bins to average> -n (calculates nyquist frequency)\n";
	 cout << "reho_csf takes a time series NIFTI volume set and runs a fourier transform over the voxels highlighted in the mask volume.\n";
	 cout << "The squared magnitude of transform result (also known as power) are averaged over a number of bins.\n";
	 cout << "Frequency domain is calculated using the total time of the input volume set.\n";
	 cout << "Averaged power results are written out to <output filename>.nii.gz as a NIFTI volume with the number of volume equal to number of bins.\n";
	 cout << "The frequency midpoint of each bin range is written to <output filename>.freq \n";
	 cout << "If bin ootions isn't used every power value is saved to the image volume set and\n\
	 it's corresponding frequency is printed instead of the bin frequency.\n";
	 cout << "Nyquist frequency option sets maximum frequency to half the sample frequency.\n";
	 cout << "Powers with a frequency greater than the nyquist frequency aren't included in output.\n";
}
int main(int argc, const char * argv []){
	if (argc == 1){
		print_help();
		return 0;
	}
	
	string input_filename ;
	string mask_filename ;
	const char * output_filename = 0;
	double time = 0.0;
	int n_bins = 0;
	bool use_nyquist = false;
	for(int i = 1 ; i < argc; i++){
		string current_argument = string(argv[i]);
		for (int j = 0 ; j < current_argument.length(); j++){
			current_argument[j] = toupper(current_argument[j]);
		}
		if((current_argument.compare("-I") == 0 || current_argument.compare("--I") == 0) && i + 1 < argc){
			input_filename = string(argv[++i]);
		}else if((current_argument.compare("-MASK") == 0 || current_argument.compare("--MASK") == 0) && i + 1 < argc){
			mask_filename = string(argv[++i]);
		}else if((current_argument.compare("-O") == 0 || current_argument.compare("--O") == 0 || \
			current_argument.compare("-OUT") == 0 || current_argument.compare("--OUT") == 0) && i + 1 < argc){
			output_filename = argv[++i];
		}else if((current_argument.compare("-T") == 0 || current_argument.compare("--T") == 0) && i + 1 < argc){
			time = atof(argv[++i]);
		}else if((current_argument.compare("-BINS") == 0 || current_argument.compare("--BINS") == 0) && i + 1 < argc){
			n_bins = atoi(argv[++i]);
		}else if((current_argument.compare("-N") == 0 || current_argument.compare("--N") == 0)){
			use_nyquist = true;			
		}else if((current_argument.compare("-HELP") == 0 || current_argument.compare("--HELP") == 0)){
			print_help();
			return 0;
		}else{
			cout << "Invalid argument was entered\n";
			return 1;
		}
	}
	/*
	if(n_bins <= 0 ){
		cout << "The number of bins either wasn't specified with -bins option or is less than equal to zero\n";
		return 1;
	}*/
	if(input_filename.length() == 0){
		cout << "No input volume was specified using the -i argument" << endl;\
		return 1;
	}

	if(mask_filename.length() == 0){
		cout << "No mask volume was specified using the -mask argument\n";
		return 1;
	}

	if(!output_filename){
		cout << "No output filename was specified using the -o argument\n";
		return 1;
	}

	if(time <= 0.0){
		cout << "Either no time was specified using argument -t or time was set to zero\n";
		return 1;
	}

	RicVolumeSet * input_volume = new RicVolumeSet(string(input_filename));
	RicVolumeSet * mask_volume = new RicVolumeSet(string(mask_filename));

	if(input_volume->nx != mask_volume->nx || input_volume->ny != mask_volume->ny || input_volume->nz != mask_volume->nz){
		cout << "The dimensions of the input volume do not match the dimensions of the mask volume\n";
		return 1;
	}
	double fs = 1.0/time;
	int n = 1;
	while (n < input_volume->nvol){
		n *= 2;
	}

	double f_domain_size = fs/n;
	
	vector<double> frequencies(n);
	int max_n = n;
	if(use_nyquist) max_n /= 2;
	for(int i = 0; i < max_n; i++){
		frequencies[i] = i*f_domain_size;
	}
	double bin_size = fs/n_bins;
	if(use_nyquist) bin_size /= 2;

	vector<double> frequency_bins;
	if(n_bins > 0){	
		frequency_bins.resize(n_bins);
		for(int i = 0 ; i < n_bins; i++){
			frequency_bins[i] = (i+1)*bin_size;
		}
	}
	string output_volume_filename = string(output_filename) + ".nii.gz";
	RicVolumeSet * output_volume;
	if(n_bins > 0){
		output_volume = new RicVolumeSet(input_volume->nx, input_volume->ny, input_volume->nz, n_bins);
	}else{
		output_volume = new RicVolumeSet(input_volume->nx, input_volume->ny, input_volume->nz, max_n);
	}
#pragma omp parallel for collapse(3) schedule(dynamic)
	for(int x = 0; x < input_volume->nx; ++x){
		for(int y = 0 ; y < input_volume->ny; ++y){
			for(int z = 0; z < input_volume->nz; ++z){
				if(mask_volume->VolSet[0].vox[x][y][z] != 0.0){
					list<complex<double>> input_list(input_volume->nvol);
					list<complex<double>>::iterator list_iter = input_list.begin();
					for(int i = 0; i < input_volume->nvol; ++i){
						*list_iter = input_volume->VolSet[i].vox[x][y][z];
						++list_iter;
					}
					input_list = pad_list(input_list);
					list<complex<double>> output_list = radix2_FFT(input_list);
					if(n_bins > 0){
						//vector<double> power_binned(n_bins);
						double power_sum = 0.0;
						unsigned power_count = 0;
						int bin_index = 0;
						list_iter = output_list.begin();
						for(int i = 0 ; i < n; ++i){

							if(frequencies[i] < frequency_bins[bin_index]){
								power_sum += norm(*list_iter++)/n;
								power_count++;
							}else{
								output_volume->VolSet[bin_index].vox[x][y][z] = power_sum/power_count;
								power_sum = norm(*list_iter++)/n;
								power_count = 1;
								bin_index++;
							}
						}
						output_volume->VolSet[bin_index].vox[x][y][z] = power_sum/power_count;
					}else{
						list_iter = output_list.begin();
						for(int i = 0 ; i < max_n; ++i){
							output_volume->VolSet[i].vox[x][y][z] = norm(*list_iter++)/n;
						}
					}
				}

			}
		}
	}
	output_volume->NIFTIorientation = input_volume->NIFTIorientation;
	
	output_volume->Write_NIFTI(output_volume_filename);
	cout << "Output NIFTI volume set written to " << output_volume_filename << endl;
	string freq_output_filename = string(output_filename) + ".freq";
	ofstream output_stream(freq_output_filename.c_str());
	if(n_bins > 0){
		output_stream << frequency_bins[0]/2.0;
		for(int i = 1; i < n_bins; ++i){
			output_stream << " " << frequency_bins[i] - bin_size/2.0;
		}
		output_stream.close();
	}else{
		output_stream << frequencies[0];
		for(int i = 1 ; i < max_n; ++i){
			output_stream << " " << frequencies[i];
		}
		output_stream.close();
	}
	cout << "Output frequency bin midpoint list written to " << freq_output_filename << endl;
	return 0;
}