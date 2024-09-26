#include <per/prCameraModelConvert.h>

//for Ceres-based optiomization (TODO:make optional)
/*
// Creates the optimization problem.
void CreateProblem(prOmni &ocam,
                   Problem* problem,
                   double *intrinsics,
                   unsigned imWidth, unsigned imHeight)
{ 
    double Xs, Ys, Zs, u, v;
    
    unsigned int nblig = 180+1; //+1 for zero deg
		double incr_phi = M_PI / 180.; // 1 degree to radian
		double phi = -180*0.5 * incr_phi;
	
    for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
		{
			Zs = cos(phi);
			Ys = 0;
			Xs = sin(phi);
      
      prPointFeature P;
      P.set_X(Xs);
      P.set_X(Ys);
      P.set_X(Zs);
            
      ocam.project3DImage(P);
      ocam.meterPixel(P);
      
      u = P.get_u();
      v = P.get_v();
      
      ceres::CostFunction* cost_function = multiSresidue::Create(Xs, Ys, Zs, u, v);
      
      problem->AddResidualBlock(cost_function,
                                NULL, // squared loss //new ceres::CauchyLoss(0.5), //
                                intrinsics,
                                );
    }
    
    //problem->SetParameterization(h, new ceres::HomogeneousVectorParameterization(9));
    
    //problem->SetParameterBlockConstant(c);
    //problem->SetParameterBlockConstant(h);
    
}

void SolveProblem(Problem* problem, double *intrinsics) {
    
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.function_tolerance = 1e-12;
    options.max_num_iterations = 500;
    options.initial_trust_region_radius = 1e1;
    
    ceres::Solver::Summary summary;
    ceres::Solve(options, problem, &summary);
    std::cout << summary.FullReport() << "\n";
    
    std::cout << "h : " << intrinsics[0] << " " << intrinsics[1] << " " << intrinsics[2] << "    " << intrinsics[3] << " " << intrinsics[4] << " " << intrinsics[5] << "    " << intrinsics[6] << " " << intrinsics[7]  << std::endl;
}


double prCameraModelConvert::convertOptBearing(prOmni &ocam, prPolyCart &pccam, unsigned imWidth, unsigned imHeight)
{
    double residual = 0;
    
    ceres::Problem problem;
    double intrinsics[7] = {ocam.getau(), ocam.getu0(), ocam.getv0(), 0., 0., 0., 0.}; //todo: make dynamic
    
    CreateProblem(ocam, problem, intrinsics, imWidth, imHeight);
    
    SolveProblem(&problem, intrinsics);
    
    double a[5] = {intrinsics[3], 0., intrinsics[4], intrinsics[5], intrinsics[6]};
    
    pccam.init(intrinsics[0], intrinsics[1], intrinsics[2], a);
    
    //residual = problem.getResidual();
    
    return residual;
}
*/

double prCameraModelConvert::convert(prOmni &ocam, prPolyCart &pccam, double theoreticalFOV, vpMatrix *abs_err)
{
	double residual = 0;
	unsigned int nblig = vpMath::round(theoreticalFOV)+1; //+1 for zero deg
	vpMatrix A(nblig, 2);
	vpColVector b(nblig);
	double incr_phi = M_PI / 180.; // 1 degree to radian
	double phi = -theoreticalFOV*0.5 * incr_phi;
	
	
	double o_au = ocam.getau();
	double o_xi = ocam.getXi();
	
	vpColVector Zs(nblig);
	vpColVector Xs(nblig);
	vpColVector up(nblig);
	
	double alpha;
	double rho_up;
	
	for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
	{
		Zs[i] = cos(phi);
		Xs[i] = sin(phi);

		up[i] = o_au*Xs[i] / (Zs[i] + o_xi);
		
		alpha += up[i] / Xs[i];
			
		//std::cout << up[i] / Xs[i] << std::endl;
	}
	
	alpha /= nblig;
	alpha = 1;
	
	for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
	{
		rho_up = sqrt(up[i]*up[i]);
	
		//A[i][0] = 1; /* A Xs[i] */ A[i][1] = rho_up*rho_up; //A[i][2] = A[i][1]*rho_up; //A[i][3] = A[i][2]*rho_up; 
		//b[i] = alpha*Zs[i];
		A[i][0] = 1; /* A Xs[i] */ A[i][1] = rho_up*rho_up; // A[i][2] = A[i][1]*rho_up; //A[i][3] = A[i][2]*rho_up; 
		b[i] = up[i]*Zs[i]/Xs[i];
	}
	
	//std::cout << A << std::endl;
	//std::cout << b << std::endl;
	
	vpColVector params = A.pseudoInverse()*b;
	
	//std::cout << params.t() << std::endl;

	double a[5] = {params[0], 0., params[1], 0., 0.};
	
	pccam.init(alpha, ocam.getu0(), ocam.getv0(), a);
	
	if(abs_err == NULL)
	{
		//compute residual: later use prCmp?
		double r_rho_up;
		for(unsigned int i = 0 ; i < nblig; i++)
		{
			//TODO: use ocam.project3DImage(...)
			r_rho_up = params[0] + params[1]*A[i][1]; //+ params[2]*A[i][2] + params[3]*A[i][3];  
			//residual += pow(Xs[i] - up[i]/alpha, 2) + pow(Zs[i] - r_rho_up/alpha, 2);
			residual += pow(r_rho_up - b[i], 2); //TODO: find another more geometric
		}
	}
	else
	{
		abs_err->resize(nblig, 2, false);
		phi = -theoreticalFOV*0.5 * incr_phi;
		//compute residual: later use prCmp?
		double r_rho_up;
		for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
		{
			(*abs_err)[i][0] = phi;
			//TODO: use ocam.project3DImage(...)
			r_rho_up = params[0] + params[1]*A[i][1]; //+ params[2]*A[i][2] + params[3]*A[i][3];  
			//residual += pow(Xs[i] - up[i]/alpha, 2) + pow(Zs[i] - r_rho_up/alpha, 2);
			(*abs_err)[i][1] = fabs(r_rho_up - b[i]);
			residual += pow((*abs_err)[i][1], 2); //TODO: find another more geometric
		}
	}
	
	residual =sqrt(residual)/(double)nblig;
	
	return residual;
}




double prCameraModelConvert::convert(prFisheyeEquidistant &fecam, prOmni &ocam, double theoreticalFOV, vpMatrix *abs_err)
{
	double residual = 0;
	unsigned int nblig = vpMath::round(theoreticalFOV)+1; //+1 for zero deg
	vpMatrix A(nblig, 2);
	vpColVector b(nblig);
	double incr_phi = M_PI / 180.; // 1 degree to radian
	double phi = -theoreticalFOV*0.5 * incr_phi;
	
	double fe_au = fecam.getau();
	
	vpColVector Zs(nblig);
	vpColVector Xs(nblig);
	vpColVector ufe(nblig);

	for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
	{
		Zs[i] = cos(phi);
		Xs[i] = sin(phi);
		ufe[i] = phi*fe_au;
	
		A[i][0] = Xs[i] / ufe[i]; A[i][1] = -1.;
		b[i] = Zs[i];
	}
	
	vpColVector params = A.pseudoInverse()*b;
	
	ocam = prOmni(params[0], params[0], fecam.getu0(), fecam.getv0(), params[1]);

	if(abs_err == NULL)
	{
		//compute residual: later use prCmp?
		for(unsigned int i = 0 ; i < nblig; i++)
		{
			//TODO: use ocam.project3DImage(...)
			residual += pow((params[0]*Xs[i]/(Zs[i]+params[1]) - ufe[i]), 2);
		}
	}
	else
	{
		abs_err->resize(nblig, 2, false);
		phi = -theoreticalFOV*0.5 * incr_phi;
		//compute residual: later use prCmp?
		for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
		{
			(*abs_err)[i][0] = phi;
			//TODO: use ocam.project3DImage(...)
			(*abs_err)[i][1] = fabs(params[0]*Xs[i]/(Zs[i]+params[1]) - ufe[i]);
			residual += pow((*abs_err)[i][1], 2);
		}
	}
	
	residual =sqrt(residual)/(double)nblig;
	
	return residual;
}

double prCameraModelConvert::convert_distortions(prCameraModel *incam, prCameraModel *outcam, double theoreticalFOV, vpMatrix *abs_err)
{
	double residual = 0;
	
	if(theoreticalFOV >= 180)
		theoreticalFOV = 90;
	

	double incr_phi = M_PI / 180.;//theoreticalFOV * 0.25; // M_PI / 180.; // 1 degree to radian
	unsigned int nblig = vpMath::round(theoreticalFOV)+1; //vpMath::round(theoreticalFOV/2)+1; //+1 for zero deg //3;//
	//theoreticalFOV *= M_PI / 180.;
	double phi = -theoreticalFOV*0.5 * incr_phi;//0.0;
	vpMatrix A(nblig, 3);
	A = 0;
	vpColVector b(nblig);
	b = 0;
	vpColVector xu(nblig);
	xu = 0;
	
	double r2, r4, r6, den;
	for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
	{
		xu[i] = sin(phi)/cos(phi);
		r2 = xu[i]*xu[i]; 
		r4 = r2*r2;
		r6 = r2*r4;

		den = 1+incam->getk6()*r2+incam->getk7()*r4+incam->getk8()*r6;
	
		A[i][0] = den*r2; A[i][1] = den*r4; A[i][2] = den*r6; 
		b[i] = 1+incam->getk1()*r2+incam->getk2()*r4+incam->getk3()*r6 - den;
	}
	
	vpColVector params = A.pseudoInverse()*b;
	
	outcam->setDistorsionParameters(params[0], params[1], params[2], incam->getk4(), incam->getk5(), 0, 0, 0);
	
	if(abs_err == NULL)
	{
		//compute residual: later use prCmp?
		for(unsigned int i = 0 ; i < nblig; i++)
		{
			//TODO: use ocam.project3DImage(...)
			prPointFeature inp, outp;
			inp.set_x(xu[i]);
			inp.set_y(0);
			incam->meterPixelConversion(inp);
			outp.set_x(xu[i]);
			outp.set_y(0);
			outcam->meterPixelConversion(outp);
			residual += pow(outp.get_u() - inp.get_u(), 2);
		}
	}
	else
	{
		abs_err->resize(nblig, 2, false);
		phi = -theoreticalFOV*0.5 * incr_phi; //0.; //

		//compute residual: later use prCmp?
		for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
		{
			//TODO: use ocam.project3DImage(...)
			prPointFeature inp, outp;
			inp.set_x(xu[i]);
			inp.set_y(0);
			incam->meterPixelConversion(inp);
			outp.set_x(xu[i]);
			outp.set_y(0);
			outcam->meterPixelConversion(outp);

			(*abs_err)[i][0] = phi;
			(*abs_err)[i][1] = fabs(outp.get_u() - inp.get_u());

			residual += pow((*abs_err)[i][1], 2);
		}
	}
	
	residual =sqrt(residual/(double)nblig);
	
	return residual;
}


