/*
 Copyright (C) 2011 Georgia Institute of Technology

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

/*
 Connor-Stevens model: Dayan and Abbott, Theoretical Neuroscience, Ch6.
 tau given in ms
 This is the Hodgkin-Huxley model plus an A-type potassium current.
 */

#include <connor-stevens.h>
#include <math.h>
#include <QtGui>

// Model Functions

static inline double
alpha_m(double V)
{
  double x = -(V + 29.7);
  double y = 10.0;

  if (fabs(x / y) < 1e-6)
    return 0.38 * y * (1 - x / y / 2.0);
  else
    return 0.38 * x / (exp(x / y) - 1.0);
}

static inline double
beta_m(double V)
{
  return 15.2 * exp(-0.0556 * (V + 54.7));
}

static inline double
m_inf(double V)
{
  return alpha_m(V) / (alpha_m(V) + beta_m(V));
}

static inline double
tau_m(double V)
{
  return 1.0e-3 / (alpha_m(V) + beta_m(V));
}

static inline double
alpha_h(double V)
{
  return 0.26 * exp(-0.04 * (V + 48.0));
}

static inline double
beta_h(double V)
{
  return 3.8 / (1.0 + exp(-0.1 * (V + 18)));
}

static inline double
h_inf(double V)
{
  return alpha_h(V) / (alpha_h(V) + beta_h(V));
}

static inline double
tau_h(double V)
{
  return 1.0e-3 / (alpha_h(V) + beta_h(V));
}

static inline double
alpha_n(double V)
{
  double x = -(V + 45.7);
  double y = 10.0;

  if (fabs(x / y) < 1e-6)
    return 0.02 * y * (1 - x / y / 2.0);
  else
    return 0.02 * x / (exp(x / y) - 1.0);
}

static inline double
beta_n(double V)
{
  return 0.25 * exp(-0.0125 * (V + 55.7));
}

static inline double
n_inf(double V)
{
  return alpha_n(V) / (alpha_n(V) + beta_n(V));
}

static inline double
tau_n(double V)
{
  return 1.0e-3 / (alpha_n(V) + beta_n(V));
}

static inline double
a_inf(double V)
{
  return pow(
      0.0761 * exp(0.0314 * (V + 94.22)) / (1 + exp(0.0346 * (V + 1.17))),
      (1 / 3.0));
}

static inline double
tau_a(double V)
{
  return 0.3632e-3 + 1.158e-3 / (1 + exp(0.0497 * (V + 55.96)));
}

static inline double
b_inf(double V)
{
  return pow(1 / (1 + exp(0.0688 * (V + 53.3))), 4.0);
}

static inline double
tau_b(double V)
{
  return 1.24e-3 + 2.678e-3 / (1 + exp(0.0624 * (V + 50.0)));
}

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
  return new ConnorStevens();
}

static DefaultGUIModel::variable_t vars[] =
  {
    { "Vm", "Membrane Potential (V)", DefaultGUIModel::OUTPUT, },
    { "Istim", "Input current (A/cm^2)", DefaultGUIModel::INPUT, },
        { "Iapp (uA/cm^2)", "Applied Current (uA/cm^2)",
            DefaultGUIModel::PARAMETER, },
        { "V0 (mV)", "Initial membrane potential (mV)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "Cm (uF/cm^2)", "Specific membrane capacitance (uF/cm^2)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "G_Na_max (mS/cm^2)", "Maximum Na+ conductance density (mS/cm^2)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "E_Na (mV)", "Sodium reversal potential (mV)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "G_K_max (mS/cm^2)",
            "Maximum delayed rectifier conductance density (mS/cm^2)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "E_K (mV)", "K+ reversal potential (mV)", DefaultGUIModel::PARAMETER
            | DefaultGUIModel::DOUBLE, },
        { "G_A_max (mS/cm^2)",
            "Maximum transient A-type K+ conductance density (mS/cm^2)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "E_A (mV)", "A-type K+ reversal potential (mV)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "G_L (mS/cm^2)", "Maximum leak conductance density mS/cm^2",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "E_L (mV)", "Leak reversal potential (mV)",
            DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
        { "Rate (Hz)", "Rate of integration (Hz)", DefaultGUIModel::PARAMETER
            | DefaultGUIModel::UINTEGER, },
        { "m", "Sodium Activation", DefaultGUIModel::STATE, },
        { "h", "Sodium Inactivation", DefaultGUIModel::STATE, },
        { "n", "Potassium Activation", DefaultGUIModel::STATE, },
        { "a", "A-type Potassium Activation", DefaultGUIModel::STATE, },
        { "b", "A-type Potassium Inactivation", DefaultGUIModel::STATE, },
        { "IKA", "A-type Potassium Current", DefaultGUIModel::STATE, },
        { "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, };

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

/*
 * Macros for making the code below a little bit cleaner.
 */

#define V (y[0])
#define m (y[1])
#define h (y[2])
#define n (y[3])
#define a (y[4])
#define b (y[5])
#define dV (dydt[0])
#define dm (dydt[1])
#define dh (dydt[2])
#define dn (dydt[3])
#define da (dydt[4])
#define db (dydt[5])
#define G_Na (G_Na_max*m*m*m*h)
#define G_K  (G_K_max*n*n*n*n)
#define G_A  (G_A_max*a*a*a*b)

ConnorStevens::ConnorStevens(void) :
  DefaultGUIModel("Connor Stevens", ::vars, ::num_vars)
{
  setWhatsThis(
      "<p><b>Connor-Stevens:</b><br>This module simulates a Connor-Stevens model neuron.</p>");
  createGUI(vars, num_vars);
  initParameters();
  update( INIT);
  refresh();
}

ConnorStevens::~ConnorStevens(void)
{
}

void
ConnorStevens::execute(void)
{
  systime = count * period; // time in seconds
  for (int i = 0; i < steps; ++i)
    solve(period / steps, y); // period in s
  output(0) = V * 1e-3; // convert to V
  count++;
}

void
ConnorStevens::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
    {
  case INIT:
    setParameter("V0 (mV)", QString::number(V0)); // initialized in mV, display in mV
    setParameter("Cm (uF/cm^2)", QString::number(Cm * 100)); // initialized in uF/mm^2, display in uF/cm^2
    setParameter("G_Na_max (mS/cm^2)", QString::number(G_Na_max * 100)); // initialized in mS/mm^2, display in mS/cm^2
    setParameter("E_Na (mV)", QString::number(E_Na)); // initialized in mV, display in mV
    setParameter("G_K_max (mS/cm^2)", QString::number(G_K_max * 100)); // initialized in mS/mm^2, display in mS/cm^2
    setParameter("E_K (mV)", QString::number(E_K)); // initialized in mV, display in mV
    setParameter("G_A_max (mS/cm^2)", QString::number(G_A_max * 100)); // initialized in mS/mm^2, display in mS/cm^2
    setParameter("E_A (mV)", QString::number(E_A)); // initialized in mV, display in mV
    setParameter("G_L (mS/cm^2)", QString::number(G_L * 100)); // initialized in mS/mm^2, display in mS/cm^2
    setParameter("E_L (mV)", QString::number(E_L)); // initialized in mV, display in mV
    setParameter("Iapp (uA/cm^2)", QString::number(Iapp * 100)); // initialized in uA/mm^2, display in uA/cm^2
    setParameter("Rate (Hz)", rate);
    setState("m", m);
    setState("h", h);
    setState("n", n);
    setState("a", a);
    setState("b", b);
    setState("IKA",IKA);
    setState("Time (s)", systime);
    break;
  case MODIFY:
    V0 = getParameter("V0 (mV)").toDouble();
    Cm = getParameter("Cm (uF/cm^2)").toDouble() / 100;
    G_Na_max = getParameter("G_Na_max (mS/cm^2)").toDouble() / 100;
    E_Na = getParameter("E_Na (mV)").toDouble();
    G_K_max = getParameter("G_K_max (mS/cm^2)").toDouble() / 100;
    E_K = getParameter("E_K (mV)").toDouble();
    G_A_max = getParameter("G_A_max (mS/cm^2)").toDouble() / 100;
    E_A = getParameter("E_A (mV)").toDouble();
    G_L = getParameter("G_L (mS/cm^2)").toDouble() / 100;
    E_L = getParameter("E_L (mV)").toDouble();
    Iapp = getParameter("Iapp (uA/cm^2)").toDouble() / 100; // displayed in uA/cm^2, calculated in uA/mm^2
    rate = getParameter("Rate (Hz)").toDouble();
    steps = static_cast<int> (ceil(period * rate));
    V = V0;
    m = m_inf(V0);
    h = h_inf(V0);
    n = n_inf(V0);
    a = a_inf(V0);
    b = b_inf(V0);
    break;
  case PERIOD:
    period = RT::System::getInstance()->getPeriod() * 1e-9; // time in seconds
    steps = static_cast<int> (ceil(period * rate));
    break;
  default:
    break;
    }
}

void
ConnorStevens::initParameters()
{
  V0 = -65; // mV
  Cm = 1e-2; // uF/mm^2
  G_Na_max = 1.2; // mS/mm^2
  G_K_max = 0.2;
  G_L = 0.003;
  G_A_max = 0.477;
  E_Na = 55.0; // mV
  E_K = -72.0;
  E_L = -70.0;
  E_A = -75.0;
  Iapp = .2404; // 1 Hz spiking
  rate = 40000;
  V = V0;
  m = m_inf(V0);
  h = h_inf(V0);
  n = n_inf(V0);
  a = a_inf(V0);
  b = b_inf(V0);
  count = 0;
  systime = 0;
  period = RT::System::getInstance()->getPeriod() * 1e-9; // s
  steps = static_cast<int> (ceil(period * rate)); // calculate how many integrations to perform per execution step
}

void
ConnorStevens::solve(double dt, double *y)
{
  double dydt[6];
  derivs(y, dydt);
  for (size_t i = 0; i < 6; ++i)
    y[i] += dt * dydt[i];
}

void
ConnorStevens::derivs(double *y, double *dydt)
{
  dV = (Iapp - input(0) * 1e6 - G_Na * (V - E_Na) - G_K * (V - E_K) - G_L * (V
      - E_L) - G_A * (V - E_A)) * 1000 / Cm;
  dm = (m_inf(V) - m) / tau_m(V);
  dh = (h_inf(V) - h) / tau_h(V);
  dn = (n_inf(V) - n) / tau_n(V);
  da = (a_inf(V) - a) / tau_a(V);
  db = (b_inf(V) - b) / tau_b(V);
  IKA = G_A * (V - E_A) * 1e-6 ; // A
}

