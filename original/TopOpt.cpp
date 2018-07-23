#include "TopOpt.h"



TopOpt::TopOpt()
{
	nelx = 20;//2
	nely = 20;//3
	volfrac = 0.4, penal = 3, rmin = 2;

	volume = volfrac * nelx*nely;

	U = MatrixXd::Zero(2 * (nely + 1)*(nelx + 1), 1);
	Ue = MatrixXd::Zero(8, 1);
	Ke = MatrixXd::Identity(8, 8);
	X = MatrixXd::Zero(nely, nelx); X.fill(volfrac);
	Xold = MatrixXd::Zero(nely, nelx);
	dc = MatrixXd::Zero(nely, nelx);;

	ke();
}

TopOpt::TopOpt(int x, int y, double v, double p, double rm)
{

	nelx = x;nely = y;
	volfrac = v, penal = p, rmin = rm;

	volume = volfrac * nelx*nely;
	U = MatrixXd::Zero(2 * (nely + 1)*(nelx + 1), 1);
	Ue = MatrixXd::Zero(8, 1);
	Ke = MatrixXd::Identity(8, 8);
	X = MatrixXd::Zero(nely, nelx); X.fill(volfrac);
	Xold = MatrixXd::Zero(nely, nelx);
	dc = MatrixXd::Zero(nely, nelx);;

	ke();
}


TopOpt::~TopOpt()
{
}

void TopOpt::ke()
{
	double E = 1; double nu = 0.3;
	Eigen::VectorXd k(9);
	k << 0.0, 0.5 - nu / 6, 0.125 + nu / 8, -0.25 - nu / 12, -0.125 + 3 * nu / 8, -0.25 + nu / 12, -0.125 - nu / 8, nu / 6, 0.125 - 3 * nu / 8;
	Ke << k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8),
		k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3),
		k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2),
		k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5),
		k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4),
		k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7),
		k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6),
		k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1);
	Ke = Ke * E / (1 - pow(nu, 2));

	//cout << Ke;


}

void TopOpt::compute()
{
	double change = 1;
	int loop = 0;
	while (change > 0.01 && loop < 1024)
	{
		cout <<"\n" <<loop << "  " << change <<"  "<<volume<< "\n";//out put per loop

		FEA();//finite element analysis

		loop += 1;
		Xold = X;
		double c = 0;//compliance
		for (int ely = 1; ely <= nely; ely++)
		{
			for (int elx = 1; elx <= nelx; elx++)
			{
				int n1 = (nely + 1)*(elx - 1) + ely;
				int n2 = (nely + 1)*elx + ely;
				Ue(1 - 1) = U(2 * n1 - 1 - 1, 1 - 1); Ue(2 - 1) = U(2 * n1 - 1, 1 - 1);
				Ue(3 - 1) = U(2 * n2 - 1 - 1, 1 - 1); Ue(4 - 1) = U(2 * n2 - 1, 1 - 1);
				Ue(5 - 1) = U(2 * n2 + 1 - 1, 1 - 1); Ue(6 - 1) = U(2 * n2 + 2 - 1, 1 - 1);
				Ue(7 - 1) = U(2 * n1 + 1 - 1, 1 - 1); Ue(8 - 1) = U(2 * n1 + 2 - 1, 1 - 1);
				double a= (Ue.transpose()*Ke*Ue)(0);
				c += pow(X(ely - 1, elx - 1), penal)*a;
				dc(ely - 1, elx - 1) = -penal * pow(X(ely - 1, elx - 1), penal - 1)*a;
			}
		}

		Filter();
		OC();
		
		change = (X - Xold).cwiseAbs().maxCoeff();
		volume = X.sum();
		
		//cout << "X\n" << X;
		//cout << "\n";
	}
	cout <<"X\n"<< X;
	cout << "\n";
	
}


void TopOpt::FEA()
{

	MatrixXd KA(2 * (nely + 1)*(nelx + 1), 2 * (nely + 1)*(nelx + 1));
	KA.setZero(); //global stiffness matrix

	MatrixXd FA(2 * (nely + 1)*(nelx + 1), 1);
	FA.setZero(); FA.coeffRef(1, 0) = -1;//set external force here
	
	for (int ely = 1; ely <= nely; ely++)
	{
		for (int elx = 1; elx <= nelx; elx++)
		{
			int n1 = (nely + 1)*(elx - 1) + ely;
			int n2 = (nely + 1)*elx + ely;

			Eigen::VectorXi num(8);
			num << 2 * n1 - 1, 2 * n1, 2 * n2 - 1, 2 * n2, 2 * n2 + 1, 2 * n2 + 2, 2 * n1 + 1, 2 * n1 + 2;
			for (int i = 1; i <= 8; i++)
			{
				for (int j = 1; j <= 8; j++)
				{
					KA.coeffRef(num(i - 1)-1, num(j - 1)-1) += pow(X(ely-1, elx-1), penal)*Ke(i - 1, j - 1);
		
				}
			}
		}
	}

	vector<int> fix; 
	fix.clear(); 
	for (int i = 1; i <= 2 * (nely + 1); i += 2)
	{
		fix.push_back(i);
	}
	fix.push_back(2 * (nelx + 1)*(nely + 1));//set fixed point

	for (int i = 1; i <= fix.size(); i++)
	{
		int ind = fix[i - 1];
		KA.col(ind - 1).setZero(); KA.row(ind - 1).setZero();
		KA.coeffRef(ind - 1, ind - 1) = 1.0;

	}	
	U= KA.llt().solve(FA);//extremely time consuming!!!
}





void TopOpt::OC()
{
	MatrixXd Xnew = X;
	double l1 = 0, l2 = 100000;
	double lmid;
	MatrixXd move = X; move.fill(0.2);
	MatrixXd temp1 = X; temp1.fill(0.001);
	while (l2 - l1 > 1e-4)
	{
		lmid = 0.5*(l1 + l2);
		MatrixXd L = X; L.fill(lmid);
		Xnew = dc.cwiseAbs().cwiseQuotient(L).cwiseSqrt().cwiseProduct(X);
		Xnew = Xnew.cwiseMin(X + move).cwiseMin(MatrixXd::Ones(nely, nelx));
		Xnew = Xnew.cwiseMax(X - move).cwiseMax(temp1);
		
		if (Xnew.sum() - volfrac * nelx*nely > 0)
		{
			l1 = lmid;
		}
		else
		{
			l2 = lmid;
		}	
	}
	X = Xnew;
}
void TopOpt::Filter()
{
	MatrixXd dcn = MatrixXd::Zero(nely, nelx);
	for (int i = 1; i <= nelx; i++)
	{
		for (int j = 1; j <= nely; j++)
		{
			double sum = 0;
			for (int k = max(int(i - rmin), 1); k <= min(int(i + rmin), nelx); k++)
			{
				for (int l = max(int(j - rmin), 1); l <= min(int(j + rmin), nely); l++)
				{
					double fac = rmin - sqrt(pow(i - k, 2) + pow(j - l, 2));
					sum += max(0.0, fac);
					dcn(j - 1, i - 1) += max(0.0, fac)*X(l - 1, k - 1)*dc(l - 1, k - 1);
				}
			}


			dcn(j - 1, i - 1) = dcn(j - 1, i - 1) / (X(j - 1, i - 1)*sum);
		}
	}
	dc = dcn;
}

void TopOpt::draw()
{
	//glClear(GL_COLOR_BUFFER_BIT);
	//glColor3f(1.0f, 0.0f, 0.0f);
	//glRectf(-0.25f, 0.25f, 0.25f, -0.25f);
	//glFlush();

	for (int i = 0; i < nely; i++)
	{
		for (int j = 0; j < nelx; j++)
		{
			float color = 1 - X(i, j);
			glClear(GL_COLOR_BUFFER_BIT);
			glColor3f(color, color, color);
			glRectf(float(-nelx + j) / 100, float(nely - i) / 100, float(-nelx + j + 1) / 100, float(nely - i - 1) / 100);
			//glFlush();
		}
	}
	//glFlush();

}