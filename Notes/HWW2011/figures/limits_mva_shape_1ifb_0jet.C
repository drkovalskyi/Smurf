{
//=========Macro generated from canvas: Canvas/Canvas
//=========  (Wed Apr 27 15:23:41 2011) by ROOT version5.28/00b
   TCanvas *Canvas = new TCanvas("Canvas", "Canvas",10,66,700,500);
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
   graph->SetPoint(0,120,3.939331);
   graph->SetPoint(1,130,2.228238);
   graph->SetPoint(2,140,1.595791);
   graph->SetPoint(3,150,1.15322);
   graph->SetPoint(4,160,0.6047384);
   graph->SetPoint(5,170,0.6825199);
   graph->SetPoint(6,180,1.037471);
   graph->SetPoint(7,190,1.713734);
   graph->SetPoint(8,200,2.139763);
   graph->SetPoint(9,250,4.55762);
   graph->SetPoint(10,300,4.606653);
   graph->SetPoint(11,300,1.085055);
   graph->SetPoint(12,250,0.8827866);
   graph->SetPoint(13,200,0.4210439);
   graph->SetPoint(14,190,0.3450514);
   graph->SetPoint(15,180,0.2055407);
   graph->SetPoint(16,170,0.1427781);
   graph->SetPoint(17,160,0.1338034);
   graph->SetPoint(18,150,0.2381324);
   graph->SetPoint(19,140,0.3138305);
   graph->SetPoint(20,130,0.491129);
   graph->SetPoint(21,120,0.8693156);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,102,318);
   Graph_Graph1->SetMinimum(0.1204231);
   Graph_Graph1->SetMaximum(5.053938);
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
   graph->SetPoint(0,120,2.621793);
   graph->SetPoint(1,130,1.547896);
   graph->SetPoint(2,140,1.136521);
   graph->SetPoint(3,150,0.7835929);
   graph->SetPoint(4,160,0.4200685);
   graph->SetPoint(5,170,0.4900784);
   graph->SetPoint(6,180,0.7096887);
   graph->SetPoint(7,190,1.146256);
   graph->SetPoint(8,200,1.482137);
   graph->SetPoint(9,250,3.128146);
   graph->SetPoint(10,300,3.27977);
   graph->SetPoint(11,300,1.548984);
   graph->SetPoint(12,250,1.364889);
   graph->SetPoint(13,200,0.640031);
   graph->SetPoint(14,190,0.4995817);
   graph->SetPoint(15,180,0.3203385);
   graph->SetPoint(16,170,0.213665);
   graph->SetPoint(17,160,0.1909262);
   graph->SetPoint(18,150,0.3606267);
   graph->SetPoint(19,140,0.5125404);
   graph->SetPoint(20,130,0.686473);
   graph->SetPoint(21,120,1.271198);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,102,318);
   Graph_Graph2->SetMinimum(0.1718336);
   Graph_Graph2->SetMaximum(3.588655);
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
   graph->SetPoint(0,120,1.955305);
   graph->SetPoint(1,130,1.123507);
   graph->SetPoint(2,140,0.8258022);
   graph->SetPoint(3,150,0.5630943);
   graph->SetPoint(4,160,0.306375);
   graph->SetPoint(5,170,0.3520466);
   graph->SetPoint(6,180,0.5128001);
   graph->SetPoint(7,190,0.8190898);
   graph->SetPoint(8,200,1.059223);
   graph->SetPoint(9,250,2.252232);
   graph->SetPoint(10,300,2.406903);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,102,318);
   Graph_Graph3->SetMinimum(0.09632223);
   Graph_Graph3->SetMaximum(2.616956);
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
   TText *text = pt->AddText("HWW #rightarrow 2l2#nu (0-jet)");
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
