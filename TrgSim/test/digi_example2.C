R__LOAD_LIBRARY(libHist)
R__LOAD_LIBRARY(libTrgSim)

void digi_example2()
{
  auto pmt = new PMTSignal();
  pmt->SetTransitTime(0);
  pmt->SetTTS(1.2);
  pmt->SetModel(PMTSignal::MOYAL);
  //pmt->SetGainScale(100.);

  //pmt->IncludeDark();
  pmt->AddHitTime(0);
  pmt->AddHitTime(30);
  pmt->Prepare();
  pmt->Draw();

  new TCanvas();
  auto fadc = new FADCWaveformGenerator();
  fadc->SetNBIT(12);
  fadc->SetVpp(200);
  fadc->SetSamplingRate(500);
  fadc->SetTimeWindow(250);
  fadc->SetPedOffset(0);
  fadc->SetPedRMS(0);

  fadc->SetSignal(pmt);
  fadc->Digitize();

  fadc->Draw();
}
