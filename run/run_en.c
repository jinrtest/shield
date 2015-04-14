#define PROTON       2212
#define NEUTRON      2112
#define APROTON     -2212
#define ANEUTRON    -2112
#define PIPLUS       211
#define PIMINUS     -211
#define PI0          111
#define KMINUS      -321
#define KPLUS        321
#define K0           130
#define AK0         -130
#define GAMMA        22
#define NUCLEI_PDG 1000000000

const Double_t PI = 4. * atan(1.);
const Double_t AngleShift = 1;
const Int_t n_hist = 13;
char names[n_hist][12] = {"Protons", "Neutrons", "A-Protons", "A-Neutrons", "Pi+", "Pi-",
                            "Pi-0", "Gamma", "K+", "K-", "K-0", "A-K-0", "Nuclei"};
const Int_t n_hist2 = 5;
char names2[n_hist2][32] = {"Nuclei", "Au (Z=79)", "Light (Z<=15)", "Medium (15<Z<=70)", "Heavy (Z>70)"};
Int_t pdgs[n_hist] = {PROTON, NEUTRON, APROTON, ANEUTRON, PIPLUS,
                      PIMINUS, PI0, GAMMA, KPLUS, KMINUS, K0, AK0, NUCLEI_PDG};
Int_t colors[n_hist] = {1, 2, 3, 4, 6, 7, 8, 12, 9, 30, 28, 46, 47};
TCanvas *c_t, *c_f;

run_en(Double_t Energy = 5, Int_t n_stat = 1000, Int_t hist_size = 200) {
    double mm = 0.1;
    cout << "root -b -q \"run_en(" << Energy << ", " << n_stat << ", " <<  hist_size << ")\"" << endl;
    cout << "::Started::" << endl;
    Double_t *hist_f_x = new Double_t[hist_size],
    **hist_f_y = new Double_t*[n_hist2],
    *hist_t_x = new Double_t[hist_size],
    **hist_t_y = new Double_t*[n_hist];

    for (Int_t j = 0; j < n_hist; j++) {
        hist_t_y[j] = new Double_t[hist_size];
    }
    for (Int_t j = 0; j < n_hist2; j++) {
        hist_f_y[j] = new Double_t[hist_size];
    }

    for (Int_t j = 0; j < hist_size; j++) {
        hist_f_x[j] = Double_t(j) / Double_t(hist_size - 1) * (PI * 0.51 - AngleShift) + AngleShift;
        for (Int_t i = 0; i < n_hist2; i++) {
            hist_f_y[i][j] = 0.0;
        }
        hist_t_x[j] = (j * PI / Double_t(hist_size - 1) - PI * 0.5) * 1.05;
        for (Int_t i = 0; i < n_hist; i++) {
            hist_t_y[i][j] = 0.0;
        }
    }
    cout << "run_en: Loading SHIELD ...";
    gSystem->Load("libGeom.so");
    gSystem->Load("libTShield.so");
    cout << " Loaded." << endl;
    cout << "run_en: Creating geometry ..." << endl;
    gROOT->ProcessLine(".x example_geometry.c");
    TShield b;

    cout << "run_en: Generation started." << endl;
    Int_t n_part;
    b.SetGeometry(gGeoManager);
    b.SetStartPoint(30*mm, 0, 0);
    for (Int_t i = 1; i <= n_stat; i++) {
        b.SetAProj(197);
        b.SetZProj(79);
        b.SetEnergy(Energy * 1e3);  // Energy in [MeV]
        Double_t angle_phi = rand() % Int_t(TMath::TwoPi() * 1E7) * 1E-7;
        Double_t angle_teta = rand() % Int_t(TMath::Pi() * 1E7) * 1E-7;
        b.SetDirection(cos(angle_teta), sin(angle_phi), cos(angle_phi));           // random should be here
        b.SetRandomSeed(clock() + time(0));
        b.GenerateEvent();

        n_part = b.GetFlyOutNum();
        cout << "Event: " << i << " / " << n_stat << " : Event generated, " << n_part << " particles imported" << endl;
        // compute angles
        TClonesArray *array = b.GetFlyOutArray();
        Int_t pdg;
        Double_t af, at;
        for (Int_t k = 0; k < n_part; k++) {
            TParticle *p = array->UncheckedAt(k);
            pdg = p->GetPdgCode();
            Double_t px = p->Px(), py = p->Py(), pz = p->Pz();
            /*if (px == 0) af = PI/2.;
            else af = atan(py/px);
            if (px < 0) af += PI;
            while (af <= 0) af += 2*PI;
            while (af > 2*PI) af -= 2*PI;*/
            if ((px == 0) && (py == 0)) {
                if (pz > 0) {
                    at = PI * 0.5;
                } else {
                    at = -PI * 0.5;
                }
            } else {
                at = atan(pz / sqrt(px * px + py * py));
            }
            for (Int_t j = 0; j < n_hist - 1; j++) { // elementary particles
                if (pdg == pdgs[j]) {
                    for (Int_t l = 0; l < hist_size - 1; l++) {
                        if ((hist_t_x[l] < at) && (hist_t_x[l + 1] >= at)) {
                            hist_t_y[j][l] += 1.0;
                        }
                    }
                }
            }
            if (pdg >= NUCLEI_PDG) { // nuclei
                for (Int_t l = 0; l < hist_size - 1; l++) {
                    if ((hist_f_x[l] < at) && (hist_f_x[l + 1] >= at)) {
                        hist_f_y[0][l] += 1.0;
                        Int_t Z = Int_t((pdg - NUCLEI_PDG) / 10000);
                        if (Z == 79) {
                            hist_f_y[1][l] += 1.0;
                        }
                        if (Z <= 15) {
                            hist_f_y[2][l] += 1.0;
                        }
                        if ((Z <= 70) && (Z > 15)) {
                            hist_f_y[3][l] += 1.0;
                        }
                        if ((Z > 70) && (Z != 79)) {
                            hist_f_y[4][l] += 1.0;
                        }
                    }
                    if ((hist_t_x[l] < at) && (hist_t_x[l + 1] >= at)) {
                        hist_t_y[n_hist - 1][l] += 1.0;
                    }
                }
            }
        }
    }

    {
        cout << "run_en: Printing dist_light_Energy" << endl;
        c_f = new TCanvas("hist_f", "Particles | angle_t (fine)", 200, 10, 600, 400);
        c_f->cd()->Divide(3, 2);
        for (Int_t i = 0; i < n_hist2; i++) {
            TGraph *gr = new TGraph(hist_size, hist_f_x, hist_f_y[i]);
            c_f->cd(i + 1)->SetLogy();
            gr->SetLineColor(colors[i]);
            gr->SetLineWidth(1);
            gr->SetLineStyle(2);
            gr->SetMarkerColor(colors[i]);
            gr->SetMarkerSize(0.6);
            gr->SetTitle(names2[i]);
            gr->Draw("Al*");
        }
        std::stringstream ss;
        ss << "dist_light_Energy" << Energy << "GeV.eps";
        c_f->Print(ss.str().c_str());
    }

    {
        cout << "run_en: Printing dist_heavy_Energy" << endl;
        c_t = new TCanvas("hist_t", "Particles | angle_t", 200, 10, 600, 400);
        c_t->cd()->Divide(5, 3);
        for (Int_t i = 0; i < n_hist; i++) {
            TGraph *gr = new TGraph(hist_size, hist_t_x, hist_t_y[i]);
            c_t->cd(i + 1)->SetLogy();
            gr->SetLineColor(colors[i]);
            gr->SetLineWidth(1);
            gr->SetLineStyle(2);
            gr->SetMarkerColor(colors[i]);
            gr->SetMarkerSize(0.6);
            gr->SetTitle(names[i]);
            gr->Draw("Al*");
        }
        std::stringstream ss;
        ss << "dist_heavy_Energy" << Energy << "GeV.eps";
        c_t->Print(ss.str().c_str());
    }

    cout << "run_en: Cleaning." << endl;
    delete hist_f_x;
    for (Int_t i = 0; i < n_hist2; i++) {
        delete hist_f_y[i];
    }
    delete hist_f_y;
    delete hist_t_x;
    for (Int_t i = 0; i < n_hist; i++) {
        delete hist_t_y[i];
    }
    delete hist_t_y;
}
