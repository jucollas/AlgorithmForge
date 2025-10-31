/*
Shamelessly taken from https://github.com/mochow13/competitive-programming-library/blob/master/Math/Fast%20Walsh-Hadamard%20Transform.cpp
and adapted by yours truly, osvarp

currently UNTESTED
*/

template <typename T>
struct FWT {
	void fwt( vector<T> & io, bool sd ) {
		for (int d = 1; d < n; d <<= 1) {
			for (int i = 0, m = d<<1; i < n; i += m) {
				for (int j = 0; j < d; j++) { /// Don't forget modulo if required
					T x = io[i+j], y = io[i+j+d];
					
					io[i+j] = (x+y), io[i+j+d] = (x-y);	// xor
					if(!sd)io[i+j]/=2,io[i+j+d]/=2;/// Modular inverse if required here
					
					// io[i+j] = (sd)?x+y:x-y;    // and
					// io[i+j+d] = (sd)?x+y:y-x; // or
				}
			}
		}
	}
	// a, b are two polynomials and n is size which is power of two
	void convolution(vector<T> &A, vector<T> &B) {
		fwt(A,1); fwt(B,1);
		for (int i = 0; i < n; i++)
			a[i] = a[i]*b[i];/// Don't forget modulo if required
		fwt(A,0);
	}
	// for a*a	
	void self_convolution(vector<T> &A) {
		fwt(A,1);
		for (int i = 0; i < n; i++)
			a[i] = a[i]*a[i];/// Don't forget modulo if required
		fwt(A,0);
	}
};
FWT<ll> fwt;