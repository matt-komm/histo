import ROOT

ROOT.gROOT.SetBatch(True)

c1 = ROOT.TCanvas("c1");

l = ROOT.TMathText()
l.SetTextAlign(23);
l.SetTextSize(0.06);
l.DrawMathText(0.50, 1.000, "$\\prod_{j\\ge0} \\left(\\sum_{k\\ge0} a_{jk}z^k\\right) = \\sum_{n\\ge0} z^n \\left(\\sum_{k_0,k_1,\\ldots\\ge0\\atop k_0+k_1+\\cdots=n} a_{0k_0}a_{1k_1} \\cdots \\right)$")
l.DrawMathText(0.50, 0.800, "W_{\\delta_1\\rho_1\\sigma_2}^{3\\beta} = U_{\\delta_1\\rho_1\\sigma_2}^{3\\beta} + {1\\over 8\\pi^2} \\int_{\\alpha_1}^{\\alpha_2} d\\alpha_2^\\prime \\left[ {U_{\\delta_1\\rho_1}^{2\\beta} - \\alpha_2^\\prime U_{\\rho_1\\sigma_2}^{1\\beta} \\over U_{\\rho_1\\sigma_2}^{0\\beta}} \\right]")

c1.Print("c1.gif")
c1.Print("c1.jpg")
c1.Print("c1.png")
c1.Print("c1.ps")
c1.Print("c1.pdf")

