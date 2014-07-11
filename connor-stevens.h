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

#include <default_gui_model.h>

/*
Connor-Stevens model: Dayan and Abbott, Theoretical Neuroscience, Ch6.
tau given in ms
This is the Hodgkin-Huxley model plus an A-type potassium current.
*/

class ConnorStevens : public DefaultGUIModel
{
	
	public:
	
		ConnorStevens(void);
		virtual~ConnorStevens(void);
	
		void execute(void);
	
	protected:
	
		void update(DefaultGUIModel::update_flags_t);
	
	private:
	
		void derivs(double *, double *);
		void solve(double, double *);
		void initParameters();
	
		double y[6];
		double period;
		int steps;
	
		double V0;
		double Cm;
		double G_Na_max;
		double E_Na;
		double G_K_max;
		double E_K;
		double G_L;
		double E_L;
		double G_A_max;
		double E_A;
		double Iapp;
		double IKA;
		double rate;
		double systime;
		long long count;
};

