#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
using namespace std;
using namespace OpenBabel;


int main(int argc, char **argv)
{

	OBSmartsPattern smarts;
	string waterPattern="[OH2]";
	smarts.Init(waterPattern);
	int c;
	char *program_name= argv[0];
	if (argc != 4) {
		cerr << " Usage: " << program_name << " <input file>  <atomic number>  <atomic number>\n";
		exit(-1);
	}
	OBMol mol;

	// Find Input filetype
	char *FileIn = NULL;
	string filename;
	//FileIn = argv[1];
	//ifstream ifs(argv[1]);
	mol.Clear();
	filename=argv[1];
	int iatom = atoi(argv[2]);
	int jatom = atoi(argv[3]);
	//OBConversion conv(&std::cin,&std::cout);
	OBConversion conv;
	OBFormat *inFormat = conv.FormatFromExt(filename.c_str());

	//if ( !conv.SetInFormat("XYZ")) {
	if (!inFormat || !conv.SetInFormat(inFormat)) {
		cerr << ": cannot read input format!" << endl;
		exit (-1);
	}
	conv.ReadFile(&mol,filename);
	if (mol.Empty()) {
		cerr << "molecule not read" << endl;
		exit(-1);
	}


	int molcount;
	double width=0.05;
	double rmax=7.00;
	int nbin=int(rmax/width);
	int* hist = new int[10000];
	double* gr = new double[10000];
	double* gs = new double[10000];


	for (int i=1;i<=10000;i++) {
		hist[i]=0;
		gr[i]=0.0;
		gs[i]=0.0;
	}

	for (;;) {
		if (mol.Empty()) 
			break;


		molcount++;
		OBAtom *atom;
		bool fragonly=false;
		vector<OBMol> mols;
		vector<OBMol>::iterator itr;
		vector<OBMol>::iterator jtr;
		OBAtom *atomi;
		OBAtom *atomj;
		mols = mol.Separate();

		double Hx,Hy,Hz,Nx,Ny,Nz;
		int iter=0;
		int numi=0;int numj=0;
		for(itr=mols.begin();itr!=mols.end();++itr) {
			bool water=smarts.Match(*itr);
			if(water) 
				continue;

			FOR_ATOMS_OF_MOL(atom,*itr)
			{
				int index=atom->GetIdx();
				//if(index!=iatom)
					//continue;

				int thisnum=atom->GetAtomicNum();
				if (thisnum!=iatom) 
					continue;

				numi++;
				Nx=atom->GetX();
				Ny=atom->GetY();
				Nz=atom->GetZ();
				//cerr << thisnum << " " << index << " " << Nx << " " << Ny << " " << Nz << endl;
				for(jtr=mols.begin();jtr!=mols.end();++jtr) {
					if(itr==jtr)
						continue;
					bool waterj=smarts.Match(*jtr);
					if(!waterj)
						continue;
					FOR_ATOMS_OF_MOL(atomj,*jtr)
					{
						int thisnum=atomj->GetAtomicNum();
						if (thisnum!=jatom) 
							continue;

						numj++;
						Hx=atomj->GetX();
						Hy=atomj->GetY();
						Hz=atomj->GetZ();

						//cerr << thisnum << " " << Hx << " " << Hy << " " << Hz << endl;
						double dist=sqrt( (Hx-Nx)*(Hx-Nx) + (Hy-Ny)*(Hy-Ny) + (Hz-Nz)*(Hz-Nz));
						iter++;
						if (dist < rmax ) {
							int bin = int(dist/width) + 1;
							//cout << bin << endl;
							hist[bin] = hist[bin] + 1;
						}
						//cout << thisidx << " " << dist << endl;
					}
					if(numj==0) {
						cerr << "atomic number " << jatom << " not found." << endl;
						exit(-1);
					}
				}
			}
			if(numi==0) {
				cerr << "atomic number " << iatom << " not found." << endl;
				exit(-1);
			}
		}
		//cerr << " --- end of frame " << molcount << " --- " << endl;

		conv.Read(&mol);
	}
	cerr << molcount << " frames collected" << endl;

	double factor = (4.0/3.0) * 3.141592653589793238 * double(molcount);
	for(int i=1;i<=nbin;i++) {
		double rupper = double(i) * width;
		double rlower = rupper - width;
		double expect = factor * (pow(rupper,3) - pow(rlower,3));
		gr[i] = double(hist[i]) / expect;
	}
	if(nbin >= 5) {
		gs[1] = (69.0*gr[1] + 4.0*gr[2] - 6.0*gr[3]
				+ 4.0*gr[4] - gr[5]) / 70.0;
		gs[2] = (2.0*gr[1] + 27.0*gr[2] + 12.0*gr[3]
				- 8.0*gr[4] + 2.0*gr[5]) / 35.0;
		for(int i = 3;i<= nbin-2;i++) {
			gs[i] = (-3.0*gr[i-2] + 12.0*gr[i-1] + 17.0*gr[i]
					+ 12.0*gr[i+1] - 3.0*gr[i+2]) / 35.0;
		}
		gs[nbin-1] = (2.0*gr[nbin-4] - 8.0*gr[nbin-3]
				+ 12.0*gr[nbin-2] + 27.0*gr[nbin-1]
				+ 2.0*gr[nbin]) / 35.0;
		gs[nbin] = (-gr[nbin-4] + 4.0*gr[nbin-3] - 6.0*gr[nbin-2]
				+ 4.0*gr[nbin-1] + 69.0*gr[nbin]) / 70.0;
		for (int i=1;i<=nbin;i++) {
			gs[i] = max(0.0,gs[i]);
		}
	}

	for (int i=1;i<=nbin;i++) {
		double dist=(double(i)-0.5)*width;
		cout << i << " " << dist << " " << gr[i] << " " << gs[i] << " " << hist[i] << endl;

	}


	return(1);

}
