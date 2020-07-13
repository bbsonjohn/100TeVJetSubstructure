{
  
   TString file1 = "events20stop.root";
   TString file2 = "events5stop.root";
   TString file3 = "events15stop.root";

   TString plot1 = "jet2MetSeparation";
   TString plot2 = "jet3MetSeparation";

   bool plot3files = true;

   TFile * f1 = TFile::Open(file1);
   if (!f1) {
      Error("hbars","file1 not found");
      return 0;
   }
   TFile * f2 = TFile::Open(file2);
   if (!f2) {
      Error("hbars","file2 not found");
      return 0;
   }
   if (plot3files) {
     TFile * f3 = TFile::Open(file3);
     if (!f3) {
        Error("hbars","file3 not found");
        return 0;
     }
   }

   TCanvas *c1 = new TCanvas("c1","c1",600,400);
   gStyle->SetOptStat(kFALSE);

   TH1F *h1a = f1->Get(plot1);
   TH1F *h1b = f2->Get(plot1);
   h1a->SetLineColor(kBlue);
   h1b->SetLineColor(kRed);
   if (plot3files) {
      TH1F *h1c = f3->Get(plot1);
      h1c->SetLineColor(kGreen);
   }

   h1a->Draw();
   h1b->Draw("same");
   if (plot3files) h1c->Draw("same");
 
   c1->Update();
   gPad->WaitPrimitive();

   TH1F *h2a = f1->Get(plot2);
   TH1F *h2b = f2->Get(plot2);
   h2a->SetLineColor(kBlue);
   h2b->SetLineColor(kRed);
   if (plot3files) {
      TH1F *h2c = f3->Get(plot2);
      h2c->SetLineColor(kGreen);
   }
   h2a->Draw();
   h2b->Draw("same");
   if (plot3files) h2c->Draw("same");
   gPad->WaitPrimitive();

   f1->Close();
   f2->Close();
   if (plot3files)  f3->Close();
   // create hint1 filled with the bins integral of h1

   delete f1;
   delete f2;

   if (plot3files) delete f3;

   // scale hint1 to the pad coordinates
/*   Float_t rightmax = 1.1*hint1->GetMaximum();
   Float_t scale = gPad->GetUymax()/rightmax;
   hint1->SetLineColor(kRed);
   hint1->Scale(scale);
   hint1->Draw("same");
   // draw an axis on the right side
   TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
   gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetLineColor(kRed);
   axis->SetTextColor(kRed);
   axis->Draw();*/
   return c1;
}

