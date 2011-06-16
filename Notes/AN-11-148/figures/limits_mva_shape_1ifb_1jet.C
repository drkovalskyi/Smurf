{
//=========Macro generated from canvas: Canvas/Canvas
//=========  (Wed Apr 27 14:47:05 2011) by ROOT version5.28/00b
   TCanvas *Canvas = new TCanvas("Canvas", "Canvas",1837,116,700,500);
   Canvas->Range(72.5,-1.25,347.5,11.25);
   Canvas->SetFillColor(0);
   Canvas->SetBorderMode(0);
   Canvas->SetBorderSize(2);
   Canvas->SetFrameBorderMode(0);
   Canvas->SetFrameBorderMode(0);
   
   TH1F *hframe__1 = new TH1F("hframe__1","",1000,100,320);
   hframe__1->SetMinimum(0);
   hframe__1->SetMaximum(10);
   hframe__1->SetDirectory(0);
   hframe__1->SetStats(0);
   hframe__1->GetXaxis()->SetTitle("Higgs mass, m_{H} [GeV/c^{2}]");
   hframe__1->GetYaxis()->SetTitle(" 95% CL Limit on #sigma/#sigma_{SM} ");
   hframe__1->Draw(" ");
   
   TGraph *graph = new TGraph(22);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ffff00");
   graph->SetFillColor(ci);

   ci = TColor::GetColor("#ffff00");
   graph->SetLineColor(ci);
   graph->SetPoint(0,120,8.309748);
   graph->SetPoint(1,130,5.943055);
   graph->SetPoint(2,140,3.23361);
   graph->SetPoint(3,150,1.835373);
   graph->SetPoint(4,160,1.099278);
   graph->SetPoint(5,170,1.331002);
   graph->SetPoint(6,180,1.838723);
   graph->SetPoint(7,190,3.149486);
   graph->SetPoint(8,200,4.155317);
   graph->SetPoint(9,250,7.601368);
   graph->SetPoint(10,300,8.005275);
   graph->SetPoint(11,300,1.763697);
   graph->SetPoint(12,250,1.460712);
   graph->SetPoint(13,200,0.8174632);
   graph->SetPoint(14,190,0.6625112);
   graph->SetPoint(15,180,0.3763135);
   graph->SetPoint(16,170,0.2948985);
   graph->SetPoint(17,160,0.2432514);
   graph->SetPoint(18,150,0.4155403);
   graph->SetPoint(19,140,0.5915743);
   graph->SetPoint(20,130,0.9176959);
   graph->SetPoint(21,120,1.779907);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,102,318);
   Graph_Graph1->SetMinimum(0.2189262);
   Graph_Graph1->SetMaximum(9.116398);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("f");
   
   graph = new TGraph(22);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   ci = TColor::GetColor("#00ff00");
   graph->SetFillColor(ci);

   ci = TColor::GetColor("#00ff00");
   graph->SetLineColor(ci);
   graph->SetPoint(0,120,5.472679);
   graph->SetPoint(1,130,3.683975);
   graph->SetPoint(2,140,2.129925);
   graph->SetPoint(3,150,1.2945);
   graph->SetPoint(4,160,0.7763572);
   graph->SetPoint(5,170,0.8869534);
   graph->SetPoint(6,180,1.28113);
   graph->SetPoint(7,190,2.147866);
   graph->SetPoint(8,200,2.83584);
   graph->SetPoint(9,250,5.240576);
   graph->SetPoint(10,300,5.574305);
   graph->SetPoint(11,300,2.480357);
   graph->SetPoint(12,250,2.367227);
   graph->SetPoint(13,200,1.261466);
   graph->SetPoint(14,190,0.9443558);
   graph->SetPoint(15,180,0.5680184);
   graph->SetPoint(16,170,0.4142176);
   graph->SetPoint(17,160,0.3468476);
   graph->SetPoint(18,150,0.6119772);
   graph->SetPoint(19,140,0.8957223);
   graph->SetPoint(20,130,1.444837);
   graph->SetPoint(21,120,2.467321);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,102,318);
   Graph_Graph2->SetMinimum(0.3121629);
   Graph_Graph2->SetMaximum(6.097051);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);
   graph->SetHistogram(Graph_Graph2);
   
   graph->Draw("f");
   TLine *line = new TLine(100,1,320,1);

   ci = TColor::GetColor("#ff0000");
   line->SetLineColor(ci);
   line->SetLineWidth(2);
   line->Draw();
   
   graph = new TGraph(11);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(2);
   graph->SetPoint(0,120,3.960197);
   graph->SetPoint(1,130,2.546984);
   graph->SetPoint(2,140,1.526075);
   graph->SetPoint(3,150,0.9465461);
   graph->SetPoint(4,160,0.5520867);
   graph->SetPoint(5,170,0.6542149);
   graph->SetPoint(6,180,0.9253502);
   graph->SetPoint(7,190,1.550574);
   graph->SetPoint(8,200,2.063635);
   graph->SetPoint(9,250,3.757173);
   graph->SetPoint(10,300,4.032746);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,102,318);
   Graph_Graph3->SetMinimum(0.2040208);
   Graph_Graph3->SetMaximum(4.380812);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);
   graph->SetHistogram(Graph_Graph3);
   
   graph->Draw("pc");
   
   TPaveText *pt = new TPaveText(0.2,0.7,0.5,0.9,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   TText *text = pt->AddText("HWW #rightarrow 2l2#nu (1-jet)");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.45,0.7,0.7,0.85,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph_Graph1","95% CL exclusion: mean","l");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph_Graph2","95% CL exclusion: 68% band","f");

   ci = TColor::GetColor("#00ff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph_Graph3","95% CL exclusion: 95% band","f");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ffff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   Canvas->Modified();
   Canvas->cd();
   Canvas->SetSelected(Canvas);
}
