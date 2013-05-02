{
  gROOT->Reset();
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L WimpSpectrum.C+");
  gROOT->ProcessLine(".x StandardVars.C");
  gStyle->SetOptStat(000000);
  
  TH2F* axes = new TH2F("axes","DarkSide-50 3 year sensitivity",
			100,5,1000,100,1e-46,1e-39);
  axes->Draw("AXIS");
  gPad->SetLogx();
  gPad->SetLogy();
  axes->GetXaxis()->SetTitle("WIMP mass [GeV]");
  axes->GetXaxis()->CenterTitle();
  axes->GetXaxis()->SetTitleOffset(1.3);
  axes->GetYaxis()->SetTitle("Cross Section [cm^{2}]");
  axes->GetYaxis()->CenterTitle();
  axes->GetYaxis()->SetTitleOffset(1.3);
  axes->SetTitle("DarkSide-50 3 year Sensitivity");
  gStyle->SetTextSize(0.04);
  //draw mssm
  TGraph* g2 = new TGraph("buchmueller.png.dat");
  g2->SetFillColor(kGray+1);
  g2->Draw("f");
  TGraph* g3 = new TGraph("buchmueller.png.dat2");
  g3->SetFillColor(kGray+2);
  g3->Draw("f");
  TText* pmssm = new TText(450,2e-45,"MSSM");
  pmssm->SetTextColor(kGray+2);
  pmssm->Draw();
  
  //draw cogent
  TGraph* gcogent = new TGraph("cogent.dat");
  gcogent->SetLineColor(kGreen);
  gcogent->SetLineWidth(2);
  gcogent->Draw("l");
  TText* pcogent = new TText(13, 3e-41,"CoGeNT");
  pcogent->SetTextColor(kGreen);
  pcogent->Draw();

  //draw dama
  TGraph* gdama1 = new TGraph("dama1.dat");
  gdama1->SetLineWidth(2);
  gdama1->SetLineColor(kMagenta);
  gdama1->Draw("l");
  TText* pdama1 = new TText(14,2e-40,"DAMA/Na");
  pdama1->SetTextColor(kMagenta);
  pdama1->Draw();
  
  TGraph* gdama2 = new TGraph("dama2.dat");
  gdama2->SetLineWidth(2);
  gdama2->SetLineColor(kMagenta);
  gdama2->Draw("l");
  TText* pdama2 = new TText(100,1e-41,"DAMA/I");
  pdama2->SetTextColor(kMagenta);
  pdama2->Draw();
  
  //draw xenon
  TGraph* gxenon1 = new TGraph("xenon.dat");
  gxenon1->SetLineColor(kBlue);
  gxenon1->SetLineWidth(2);
  gxenon1->SetLineStyle(7);
  gxenon1->Draw("l");
  TText* pxenon = new TText(400,1.8e-44,"Xenon100");
  pxenon->SetTextColor(kBlue);
  pxenon->Draw();
  
  //draw cdms
  TGraph* gcdms1 = new TGraph("cdms.dat");
  gcdms1->SetLineColor(kRed);
  gcdms1->SetLineWidth(2);
  gcdms1->SetLineStyle(9);
  gcdms1->Draw("l");
  TText* pcdms = new TText(400,3e-43,"CDMSII");
  pcdms->SetTextColor(kRed);
  pcdms->Draw();
  
  TRolke rolke;
  rolke.SetPoissonBkgKnownEff(0,1,20,0.5);
  double counts = rolke.GetUpperLimit();
  counts = 4.6;
  double masstime = 100*kg*year;
  
  //double counts = 2.3;
  rho = 0.3*GeV/cm3;
  v0 = 220*km/s;
  ve = 244*km/s;
  //vesc = 544*km/s;

  double emin = 30, emax=200;
  double A=AAr;
  
  
  cout<<"The 90% upper CL on mean counts is "<<counts<<endl;
  TF1* f1 = new TF1("f1","DifferentialWimpRate(x*[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])",0,200);
  
  f1->SetParameters(keV,sigma,rho,A,18,1,Md,vesc,v0,ve);
  
  vector<double> x,y;
  for(double mass = 1; mass<10000; mass *=1.01){
    f1->SetParameter(6,mass*GeV);
    double rate = masstime*f1->Integral(emin,emax);
    if(rate==0) continue;
    double limit = sigma*counts/rate;
    //cout<<mass<<" "<<rate<<" "<<limit<<endl;
    x.push_back(mass);
    y.push_back(limit);
  }
  
  TGraph* g = new TGraph(x.size(), &x[0], &y[0]);
  g->Draw("l");
  g->SetLineWidth(2);
  TText* t = new TText(22,1.2e-45,"DarkSide-50");
  t->Draw();
  
  gPad->Update();
   
  cout<<"Sensitivity to 50 GeV WIMP is "<<g->Eval(50)<<" cm2"<<endl;
  cout<<"Sensitivity to 100 GeV WIMP is "<<g->Eval(100)<<" cm2"<<endl;
  
  if(0){
    //draw xenon
    rolke.SetGaussBkgGaussEff(3,1.8,0.32, 0.03*0.32, 0.6);
    counts = rolke.GetUpperLimit();
    masstime = 4848*kg*day;
    emin = 8.4; 
    emax=44.6;
    A=AXe;  
    f1->SetParameters(keV,sigma,rho,A,1,1,Md,vesc,v0,ve);
    x.clear(); y.clear();
    for(double mass = 1; mass<10000; mass *=1.01){
      f1->SetParameter(6,mass*GeV);
      double rate = masstime*f1->Integral(emin,emax);
      if(rate==0) continue;
      double limit = sigma*counts/rate;
      //cout<<mass<<" "<<rate<<" "<<limit<<endl;
      x.push_back(mass);
      y.push_back(limit);
    }
    TGraph* gxe = new TGraph(x.size(), &x[0], &y[0]);
    gxe->SetName("gxe");
    gxe->SetTitle("xenon");
    gxe->SetLineWidth(2);
    gxe->SetLineStyle(7);
    gxe->SetLineColor(kBlue);
    gxe->Draw("l");
    
    //draw cdms
    rolke.SetGaussBkgGaussEff(2,0.9, 0.32, 0.02, 0.2);
    counts = rolke.GetUpperLimit();
    masstime = 612*kg*day;
    //how to get it to scale properly???
    masstime *= 2.5;  
    emin = 10; 
    emax=100;
    A=AGe;  
    f1->SetParameters(keV,sigma,rho,A,18,1,Md,vesc,v0,ve);
    x.clear(); y.clear();
    for(double mass = 1; mass<10000; mass *=1.01){
      f1->SetParameter(6,mass*GeV);
      double rate = masstime*f1->Integral(emin,emax);
      if(rate==0) continue;
      double limit = sigma*counts/rate;
      //cout<<mass<<" "<<rate<<" "<<limit<<endl;
      x.push_back(mass);
      y.push_back(limit);
    }
    TGraph* gcdms = new TGraph(x.size(), &x[0], &y[0]);
    gcdms->SetName("gcdms");
    gcdms->SetTitle("cdms");
    gcdms->SetLineWidth(2);
    gcdms->SetLineStyle(9);
    gcdms->SetLineColor(kRed);
    gcdms->Draw("l");
  }
  
 


}
