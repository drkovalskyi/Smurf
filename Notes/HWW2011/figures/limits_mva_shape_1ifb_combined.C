{
//=========Macro generated from canvas: Canvas/Canvas
//=========  (Wed Apr 27 22:58:38 2011) by ROOT version5.28/00b
   TCanvas *Canvas = new TCanvas("Canvas", "Canvas",10,32,700,500);
   Canvas->Range(72.5,-0.625,347.5,5.625);
   Canvas->SetFillColor(0);
   Canvas->SetBorderMode(0);
   Canvas->SetBorderSize(2);
   Canvas->SetFrameBorderMode(0);
   Canvas->SetFrameBorderMode(0);
   
   TH1F *hframe__1 = new TH1F("hframe__1","",1000,100,320);
   hframe__1->SetMinimum(0);
   hframe__1->SetMaximum(5);
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
   graph->SetPoint(0,120,4.419194);
   graph->SetPoint(1,130,3.067068);
   graph->SetPoint(2,140,1.971036);
   graph->SetPoint(3,150,1.305519);
   graph->SetPoint(4,160,0.7214675);
   graph->SetPoint(5,170,0.6510104);
   graph->SetPoint(6,180,0.8818447);
   graph->SetPoint(7,190,2.054659);
   graph->SetPoint(8,200,2.728464);
   graph->SetPoint(9,250,5.012154);
   graph->SetPoint(10,300,5.851736);
   graph->SetPoint(11,300,0.5261396);
   graph->SetPoint(12,250,0.3345498);
   graph->SetPoint(13,200,0.1447748);
   graph->SetPoint(14,190,0.131539);
   graph->SetPoint(15,180,0.1227365);
   graph->SetPoint(16,170,0.1065977);
   graph->SetPoint(17,160,0.06793266);
   graph->SetPoint(18,150,0.1145255);
   graph->SetPoint(19,140,0.1271123);
   graph->SetPoint(20,130,0.2789323);
   graph->SetPoint(21,120,0.5683139);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,102,318);
   Graph_Graph1->SetMinimum(0.06113939);
   Graph_Graph1->SetMaximum(6.430116);
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
   graph->SetPoint(0,120,2.544547);
   graph->SetPoint(1,130,1.428868);
   graph->SetPoint(2,140,1.0466);
   graph->SetPoint(3,150,0.6873797);
   graph->SetPoint(4,160,0.3825375);
   graph->SetPoint(5,170,0.4181538);
   graph->SetPoint(6,180,0.5767955);
   graph->SetPoint(7,190,1.059875);
   graph->SetPoint(8,200,1.369473);
   graph->SetPoint(9,250,2.588085);
   graph->SetPoint(10,300,2.899472);
   graph->SetPoint(11,300,0.8984383);
   graph->SetPoint(12,250,0.7481628);
   graph->SetPoint(13,200,0.3373718);
   graph->SetPoint(14,190,0.2779389);
   graph->SetPoint(15,180,0.2368189);
   graph->SetPoint(16,170,0.1694231);
   graph->SetPoint(17,160,0.1242643);
   graph->SetPoint(18,150,0.2276617);
   graph->SetPoint(19,140,0.2831112);
   graph->SetPoint(20,130,0.4972845);
   graph->SetPoint(21,120,1.016029);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,102,318);
   Graph_Graph2->SetMinimum(0.1118379);
   Graph_Graph2->SetMaximum(3.176993);
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
   graph->SetPoint(0,120,1.829753);
   graph->SetPoint(1,130,1.020349);
   graph->SetPoint(2,140,0.678648);
   graph->SetPoint(3,150,0.4643194);
   graph->SetPoint(4,160,0.2661571);
   graph->SetPoint(5,170,0.2976494);
   graph->SetPoint(6,180,0.4120619);
   graph->SetPoint(7,190,0.678622);
   graph->SetPoint(8,200,0.8722544);
   graph->SetPoint(9,250,1.699419);
   graph->SetPoint(10,300,1.944841);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,102,318);
   Graph_Graph3->SetMinimum(0.09828861);
   Graph_Graph3->SetMaximum(2.11271);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);
   graph->SetHistogram(Graph_Graph3);
   
   graph->Draw("lp");
   
   TPaveText *pt = new TPaveText(0.2,0.7,0.5,0.9,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   TText *text = pt->AddText("HWW #rightarrow 2l2#nu");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.45,0.7,0.7,0.85,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph","95% CL exclusion: mean","l");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph","95% CL exclusion: 68% band","f");

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
   entry=leg->AddEntry("Graph","95% CL exclusion: 95% band","f");

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
