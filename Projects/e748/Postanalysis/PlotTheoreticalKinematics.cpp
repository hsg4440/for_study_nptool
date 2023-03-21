#include "NPReaction.h"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TString.h>

void PlotTheoreticalKinematics()
{
    double beam10Be {28. * 10.013534};
    double beam12Be {30. * 12.026922};
    //simply compare reaction kinematics
    NPL::Reaction r10Be(TString::Format("10Be(d,3He)9Li@%.f", beam10Be).Data());
    NPL::Reaction r12Be (TString::Format("12Be(d,3He)11Li@%.f", beam12Be).Data());

    //get graphs
    auto* g9Li {r10Be.GetKinematicLine3()};
    auto* g11Li {r12Be.GetKinematicLine3()};
    //also Elab vs thetaCM
    auto* gcm9Li {r10Be.GetELabVersusThetaCM()};
    auto* gcm11Li {r12Be.GetELabVersusThetaCM()};
    //theta3 vs theta4 in LAB
    auto* gtheta9Li {r10Be.GetTheta3VsTheta4()};
    auto* gtheta11Li {r12Be.GetTheta3VsTheta4()};
    
    //plotting
    auto* c1 {new TCanvas("c1", "Kinematic lines comparaison")};
    c1->DivideSquare(4);
    c1->cd(1);
    g11Li->SetTitle(";#theta_{3}^{LAB} [degree];T_{3} [MeV]");
    g9Li->SetLineWidth(2); g9Li->SetLineColor(kRed);
    g11Li->SetLineWidth(2); g11Li->SetLineColor(kBlue);
    g11Li->Draw("apl");
    g9Li->Draw("pl same");

    c1->cd(2);
    gcm11Li->SetTitle(";#theta_{3}^{CM} [degree];T_{3} [MeV]");
    gcm9Li->SetLineWidth(2); gcm9Li->SetLineColor(kRed);
    gcm11Li->SetLineWidth(2); gcm11Li->SetLineColor(kBlue);
    gcm11Li->Draw("apl");
    gcm9Li->Draw("pl same");

    c1->cd(3);
    gtheta9Li->SetTitle(";#theta_{3}^{LAB} [degree];#theta_{4}^{LAB} [degree]");
    gtheta9Li->SetLineWidth(2); gtheta9Li->SetLineColor(kRed);
    gtheta11Li->SetLineWidth(2); gtheta11Li->SetLineColor(kBlue);
    gtheta9Li->Draw("apl");
    gtheta11Li->Draw("pl same");
}
