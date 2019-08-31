#include <objective_functions/total_objective.h>

TotalObjective::TotalObjective()
{

}

void TotalObjective::init()
{
	//assume that each objective's member have been set outside
	for (auto &objective : objectiveList) {
		objective->init();
	}
	prepare_hessian();
}

void TotalObjective::updateX(const VectorXd& X)
{
	for (auto &objective : objectiveList)
		objective->updateX(X);
}

double TotalObjective::value()
{
	double f=0;
    for (auto &objective : objectiveList)
        if (objective->w != 0)
		    f+= objective->w*objective->value();

	return f;
}

void TotalObjective::gradient(VectorXd& g)
{
	VectorXd gi;
	g.setZero();
	for (auto &objective : objectiveList) {
		if (objective->w != 0)
		{
			objective->gradient(gi);
			g += objective->w*gi;
		}
	}
}

void TotalObjective::hessian()
{
	SS.clear();
	for (auto const &objective : objectiveList)
	{
        if (objective->w != 0) //Just don't update the hessian, but we still must enter those elements into the big hessian to have the same sparsity pattern
		    objective->hessian();
        vector<double> SSi; SSi.resize(objective->SS.size());
        for (int i = 0; i < objective->SS.size(); i++)
            SSi[i] = objective->w * objective->SS[i];
		SS.insert(SS.end(), SSi.begin(), SSi.end());
	}
}

void TotalObjective::prepare_hessian()
{
	//assume that each subobjective already prepared its hessian
	II.clear(); JJ.clear(); SS.clear();
	for (auto &objective : objectiveList)
	{
		II.insert(II.end(), objective->II.begin(), objective->II.end());
		JJ.insert(JJ.end(), objective->JJ.begin(), objective->JJ.end());
		SS.insert(SS.end(), objective->SS.begin(), objective->SS.end());
	}
}