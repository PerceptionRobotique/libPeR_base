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
//UC2CP
double prCameraModelConvert::convert(prOmni &ocam, prPolyCart &pccam, double theoreticalFOV, vpMatrix *abs_err)
{
    double residual = 0;
    unsigned int nblig = vpMath::round(theoreticalFOV) + 1; // +1 for zero deg
    vpMatrix A(nblig, 2);
    vpColVector b(nblig);
    double incr_phi = M_PI / 180.; // 1 degree to radian
    double phi = -theoreticalFOV * 0.5 * incr_phi;

    double o_au = ocam.getau();
    double o_xi = ocam.getXi();

    vpColVector Zs(nblig);
    vpColVector Xs(nblig);
    vpColVector up(nblig);

    double alpha;
    double rho_up;

    for (unsigned int i = 0; i < nblig; i++, phi += incr_phi)
    {
        Zs[i] = std::cos(phi);
        Xs[i] = sin(phi);

        up[i] = o_au * Xs[i] / (Zs[i] + o_xi);

        alpha += up[i] / Xs[i];
    }

    alpha /= nblig;
    alpha = 1;

    for (unsigned int i = 0; i < nblig; i++, phi += incr_phi)
    {
        rho_up = sqrt(up[i] * up[i]);

        A[i][0] = 1;
        A[i][1] = rho_up * rho_up;
        b[i] = up[i] * Zs[i] / Xs[i];
    }

    vpColVector params = A.pseudoInverse() * b;
    double a[5] = {params[0], 0., params[1], 0., 0.};

	pccam.init(alpha, ocam.getu0(), ocam.getv0(), a);
	double u_inp,x_inp;	
	double u0 = ocam.getu0();	
    if (abs_err == NULL)
    {
        double r_rho_up;
        for (unsigned int i = 0; i < nblig; i++)
        {
            // Create a prPointFeature for the 3D point (Xs[i], Zs[i])
            prPointFeature inp, outp;
            inp.set_X(Xs[i]);
            inp.set_Y(0);
            inp.set_Z(Zs[i]);
			ocam.project3DImage(inp);

			x_inp = inp.get_x();
            u_inp = alpha*x_inp + u0;
			
			r_rho_up = inp.get_x(); 
            residual += pow(r_rho_up , 2); 
        }
    }
    else
	{
		abs_err->resize(nblig, 2, false);
		phi = -theoreticalFOV * 0.5 * incr_phi;
		double phi_out, X_out, Z_out,phi_in,X_inp, u,x, x_out, u_out ;
		
		for (unsigned int i = 0; i < nblig; i++, phi += incr_phi)
		{
			(*abs_err)[i][0] = phi;
			prPointFeature inp, outp;

			//input
			inp.set_X(Xs[i]);
			inp.set_Y(0);
			inp.set_Z(Zs[i]);
			//ocam.project3DImage(inp);
			x = Xs[i]/(Zs[i] + ocam.getXi());
			u = o_au * x + u0;
			phi_in = phi * 180./(M_PI); 

			//output
			X_out = u - u0;
			Z_out = params[0] + params[1] * std::pow((u - u0),2);
						
			x_out = X_out/(Z_out + ocam.getXi()*std::sqrt(X_out*X_out+Z_out*Z_out));
			u_out = o_au * x_out + u0;

			phi_out = std::atan2(X_out, Z_out) * 180./(M_PI);

			(*abs_err)[i][0] = phi;
			(*abs_err)[i][1] = std::fabs(u-u_out);

			residual += std::pow((*abs_err)[i][1], 2);
		}
	}


    residual = sqrt(residual / (double)nblig);

    return residual;
}

//EF2UC
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
//RP2P
double prCameraModelConvert::convert_distortions(prCameraModel *incam, prCameraModel *outcam, double theoreticalFOV, vpMatrix *abs_err)// vpMatrix *abs_err)
{
	double residual = 0;
	
	if(theoreticalFOV >= 180)
		theoreticalFOV = 90;
	

	double incr_phi = M_PI / 180.;//theoreticalFOV * 0.25; // M_PI / 180.; // 1 degree to radian
	unsigned int nblig = vpMath::round(theoreticalFOV)+1; //vpMath::round(theoreticalFOV/2)+1; //+1 for zero deg //3;//
	//theoreticalFOV *= M_PI / 180.;
	double phi = -theoreticalFOV*0.5 * incr_phi;;
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
//RP2UC
double prCameraModelConvert::convert(prCameraModel *incam, prOmni *outcam, double theoreticalFOV, vpMatrix *abs_err)
{
	double residual = 0;

	double incr_phi = M_PI / 180.;
	unsigned int nblig = vpMath::round(theoreticalFOV)+1; 
	double phi = -theoreticalFOV*0.5 * incr_phi;
	vpMatrix A(nblig, 2);
	A = 0;
	vpColVector b(nblig);
	b = 0;
	vpColVector xu(nblig);
	xu = 0;
	
	std::vector<prPointFeature> v_pt;
	double r2, r4, r6, den;
	for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
	{
		double Zs = cos(phi);
		double Xs = sin(phi);
		prPointFeature pt;
		pt.set_X(Xs);
		pt.set_Y(0);
		pt.set_Z(Zs);
		incam->project3DImage(pt);
		incam->meterPixelConversion(pt);
		v_pt.push_back(pt);
		A[i][0] = pt.get_u() - incam->getu0();//sin(phi)/cos(phi);
		A[i][1] = -Xs;
	
		b[i] = -Zs*A[i][0];
	}
	
	vpColVector params = A.pseudoInverse()*b;
	
	outcam->init(params[1], params[1], incam->getu0(), incam->getv0(), params[0]);
	
	if(abs_err == NULL)
	{
		//compute residual: later use prCmp?
		for(unsigned int i = 0 ; i < nblig; i++)
		{
			//TODO: use ocam.project3DImage(...)
			prPointFeature outp;
			outp = v_pt[i];
			outcam->project3DImage(outp);
			outcam->meterPixelConversion(outp);
			residual += pow(outp.get_u() - v_pt[i].get_u(), 2);
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
			prPointFeature outp;
			outp = v_pt[i];
			outcam->project3DImage(outp);
			outcam->meterPixelConversion(outp);

			(*abs_err)[i][0] = phi;
			(*abs_err)[i][1] = fabs(outp.get_u() - v_pt[i].get_u());

			residual += pow((*abs_err)[i][1], 2);
		}
	}
	
	residual = sqrt(residual/(double)nblig);
	
	return residual;
}

//RP2OF
double prCameraModelConvert::conversion(prCameraModel *inputcam, prCameraModel *outputcam, double theoreticalFOV, vpMatrix *abs_err) 
{
    double residual = 0;
    
    if (theoreticalFOV >= 180) {
        theoreticalFOV = 90;
    }

    double incr_phi = M_PI / 180.; // 1 degree to radian
    unsigned int nblig = vpMath::round(theoreticalFOV) + 1;

    vpColVector Z(nblig);
    vpColVector Y(nblig);
    vpColVector X(nblig);

    double phi = -theoreticalFOV * 0.5 * incr_phi;

    // Sweep of X and Z with Y set to 0. Rotation around Z 
    double rot = std::sqrt(2.0) / 2.0;

    for (unsigned int i = 0; i < nblig; i++, phi += incr_phi) {
        X[i] = rot*sin(phi);
        Y[i] = X[i];
        Z[i] = cos(phi);
    }

    vpMatrix A(nblig * 2, 6);
    vpColVector b(nblig * 2);

    for (unsigned int i = 0, j = 0; i < nblig; i++, j += 2) {
 
        double x = X[i] / Z[i];
        double y = Y[i] / Z[i];
		double rho_o = sqrt(x*x+y*y);
		double rho_x = sqrt(x*x+y*y) ;
        double theta_i 	= std::atan(rho_o);
        double theta3_i = std::pow(theta_i, 3);
        double theta5_i = std::pow(theta_i, 5);
        double theta7_i = std::pow(theta_i, 7);
        double theta9_i = std::pow(theta_i, 9);

		double alpha_u = inputcam->getau(); 
		double alpha_v = inputcam->getav(); 

		//Distortions
        double d_R = (1 + inputcam->getk1() * std::pow(rho_x, 2) + inputcam->getk2() * std::pow(rho_x, 4) + inputcam->getk3() * std::pow(rho_x, 6)) /
                     (1 + inputcam->getk6() * std::pow(rho_x, 2) + inputcam->getk7() * std::pow(rho_x, 4) + inputcam->getk8() * std::pow(rho_x, 6));

        double t_x = 2 * inputcam->getk4() * x * y + inputcam->getk5() * (std::pow(rho_x, 2) + 2 * std::pow(x, 2));
        double t_y = 2 * inputcam->getk5() * x * y + inputcam->getk4() * (std::pow(rho_x, 2) + 2 * std::pow(y, 2));

		double term1 = -alpha_u * rho_o * (d_R * x + t_x) / x;
		double term2 = -alpha_v * rho_o * (d_R * y + t_y) / y;

		//matrix
        A[j][0] = term1;
        A[j][1] = 0;
        A[j][2] = theta3_i;
        A[j][3] = theta5_i;
        A[j][4] = theta7_i;
        A[j][5] = theta9_i;

        A[j+1][0] = 0;
        A[j+1][1] = term2;
		A[j+1][2] = theta3_i;
        A[j+1][3] = theta5_i;
        A[j+1][4] = theta7_i;
        A[j+1][5] = theta9_i;

        b[j] = -theta_i;
        b[j+1] = -theta_i;
    }

    // Solve for the parameters using pseudo-inverse
    vpColVector params = A.pseudoInverse() * b;

    // Set the distortion parameters in the output camera model
	double f_x = 1/params[0];
	double f_y = 1/params[1];
	outputcam->setDistorsionParameters(params[2], params[3], params[4], params[5], 0 , 0, 0, 0);
	std::cout << "f_x: " << f_x << " and f_y: " << f_y << std::endl;

    // Compute and return the residual (if necessary)
    if (abs_err == NULL) {
        for (unsigned int i = 0; i < nblig; i++) {
            prPointFeature inp, outp;
            inp.set_x(X[i]/Z[i]);
            inp.set_y(Y[i]/Z[i]);
            inputcam->meterPixelConversion(inp);
            outp.set_x(X[i]/Z[i]);
            outp.set_y(Y[i]/Z[i]);
			double x = outp.get_x(), y = outp.get_y(), r, xd, yd;
			double u0 = inputcam->getu0(), v0 = inputcam->getv0();
			double rho_o = sqrt(x*x+y*y);
        	double theta_i 	= std::atan(rho_o);
			double theta3_i = std::pow(theta_i, 3);
			double theta5_i = std::pow(theta_i, 5);
			double theta7_i = std::pow(theta_i, 7);
			double theta9_i = std::pow(theta_i, 9);

			r = theta_i + params[2]*theta3_i + params[3]*theta5_i + params[4]*theta7_i + params[5]*theta9_i;

			xd = r * x / rho_o;
			yd = r * y / rho_o;

			outp.set_u(f_x * xd + u0);
			outp.set_v(f_y * yd + v0);

            residual += std::pow(outp.get_v() - inp.get_v(), 2);
        }
    } else {
        abs_err->resize(nblig, 2, false);
        phi = -theoreticalFOV * 0.5 * incr_phi;

        for (unsigned int i = 0; i < nblig; i++, phi += incr_phi) {
            prPointFeature inp, outp;
            inp.set_x(X[i]/Z[i]);
            inp.set_y(Y[i]/Z[i]);
            inputcam->meterPixelConversion(inp);
            outp.set_x(X[i]/Z[i]);
            outp.set_y(Y[i]/Z[i]);
			double x = outp.get_x(), y = outp.get_y(), r, xd, yd;
			double u0 = inputcam->getu0(), v0 = inputcam->getv0();
			double rho_o = sqrt(x*x+y*y);
        	double theta_i 	= std::atan(rho_o);
			double theta3_i = std::pow(theta_i, 3);
			double theta5_i = std::pow(theta_i, 5);
			double theta7_i = std::pow(theta_i, 7);
			double theta9_i = std::pow(theta_i, 9);

			r = theta_i + params[2]*theta3_i + params[3]*theta5_i + params[4]*theta7_i + params[5]*theta9_i;

			xd = r * x /rho_o;
			yd = r * y / rho_o;

			outp.set_u(f_x * xd + u0);
			outp.set_v(f_y * yd + v0);

            (*abs_err)[i][0] = phi;
            (*abs_err)[i][1] = std::fabs(outp.get_v() - inp.get_v());

            residual += std::pow((*abs_err)[i][1], 2);
        }
    }

    residual = sqrt(residual / (double)nblig);

    return residual;
}

//CP2OF
double prCameraModelConvert::FCVreverseConvert(prPolyCart *input, prCameraModel *output, double theoreticalFOV, vpMatrix *abs_err, double image_size) 
{
	double residual = 0;
	double nblig = image_size; 
	vpColVector x(nblig);
	vpColVector Up(nblig);

	double up = - nblig/2; //u-u0
	double a_0 = input->geta0();
	double a_2 = input->geta2();
	double a_3 = input->geta3();
	double a_4 = input->geta4();


	double sign_chg = -1; 
	bool is_negative_start = false; 

	for (unsigned int i = 0; i < nblig; i++, up++) {
		if (up == 0) {up++;}

		double rho_up_2 = std::pow(std::abs(up), 2);
		double rho_up_3 = std::pow(std::abs(up), 3);
		double rho_up_4 = std::pow(std::abs(up), 4);
		x[i] = up  / (a_0 + a_2 * rho_up_2 + a_3 * rho_up_3 + a_4 * rho_up_4 );

		if (i == 0){
			//std::cout << "x[0] = " << x[i] << std::endl;
			if (x[i] < 0) {
				is_negative_start = true;
				}
			//std::cout << std::boolalpha << "is_negative_start set to " << is_negative_start << std::endl;
		}

		if (i > 0 && (x[i] * x[i-1] < 0)) { // Change of sign
			if (sign_chg == -1) { 
				sign_chg = i;
			}
		}
		Up[i] = up;
	}

	// If we start with negative values, nblig_crop takes the value of nblig
	double nblig_crop = x.size()-2*sign_chg;

	if (is_negative_start) {
		nblig_crop = nblig;
		//std::cout << "And nblig_crop set equal to nblig = " << nblig << std::endl;
		sign_chg = 0;
	}

	vpMatrix A(nblig_crop, 5);
	vpColVector b(nblig_crop);
	for (unsigned int i = 0; i < nblig_crop; i++) {
		double offset = sign_chg+i;
		double theta_i 	= std::atan(std::abs(x[offset]));
		double theta3_i = std::pow(theta_i, 3);
		double theta5_i = std::pow(theta_i, 5);		
		double theta7_i = std::pow(theta_i, 7);
		double theta9_i = std::pow(theta_i, 9);
		
		A[i][0] = -Up[offset] * x[offset] / std::abs(x[offset]);
		A[i][1] = theta3_i;
		A[i][2] = theta5_i;		
		A[i][3] = theta7_i;
		A[i][4] = theta9_i;

		b[i] = -theta_i; //( Up[offset] * std::abs(x[offset]) / x[offset] );
	}
	
	vpColVector params = A.pseudoInverse() * b;
	double f_x = 1/params[0];
	double k_o1 = params[1];
	double k_o2 = params[2];
	double k_o3 = params[3];
	double k_o4 = params[4];
	
	output->setDistorsionParameters(k_o1, k_o2, k_o3, k_o4);
	output->setPixelRatio(f_x,f_x);

	if (abs_err == NULL) {
		for (unsigned int i = sign_chg; i < nblig-sign_chg; i++) {
			double offset = i;
			prPointFeature inp, outp;
			inp.set_u(Up[offset]);
			double u0 = input->getu0(),  f_x = params[0];

			double Xs,Zs;
			double up = Up[offset];
			double rho_up_2 = std::pow(std::abs(up), 2);
			double rho_up_3 = std::pow(std::abs(up), 3);
			double rho_up_4 = std::pow(std::abs(up), 4);
			Xs = up;
			Zs = a_0 + a_2 * rho_up_2 + a_3 * rho_up_3 + a_4 * rho_up_4;
			//input->projectImageSphere(inp,Xs,Ys,Zs);

			double x_FCV = Xs/Zs;
			double r, xd;
			double rho_o = std::abs(x_FCV);
			double theta_i 	= std::atan(rho_o);
			double theta3_i = std::pow(theta_i, 3);
			double theta5_i = std::pow(theta_i, 5);
			double theta7_i = std::pow(theta_i, 7);
			double theta9_i = std::pow(theta_i, 9);

			r = theta_i + k_o1*theta3_i + k_o2*theta5_i + k_o3*theta7_i + k_o4*theta9_i;

			xd = r * x_FCV / rho_o;

			outp.set_u(f_x * xd );

			residual += std::pow(outp.get_u() - inp.get_u(), 2);
		}

	} else {
		double count = 0;
		abs_err->resize(nblig_crop, 2, false);
		for (unsigned int i = 0; i < nblig_crop; i++) {
			count++;
			double offset = sign_chg + i;
			prPointFeature inp, outp;
			inp.set_u(Up[offset]);
			double Xs,Zs;
			double up = Up[offset];
			double rho_up_2 = std::pow(std::abs(up), 2);
			double rho_up_3 = std::pow(std::abs(up), 3);
			double rho_up_4 = std::pow(std::abs(up), 4);
			Xs = up;
			Zs = a_0 + a_2 * rho_up_2 + a_3 * rho_up_3 + a_4 * rho_up_4;
			//input->projectImageSphere(inp,Xs,Ys,Zs);

			double x_FCV = Xs/Zs;
			double r, xd;
			double rho_o = std::abs(x_FCV);
			double theta_i 	= std::atan(std::abs(x_FCV));
			double theta3_i = std::pow(theta_i, 3);
			double theta5_i = std::pow(theta_i, 5);
			double theta7_i = std::pow(theta_i, 7);
			double theta9_i = std::pow(theta_i, 9);

			r = theta_i + k_o1*theta3_i + k_o2*theta5_i + k_o3*theta7_i + k_o4*theta9_i;

			xd = r * x_FCV / rho_o;

			outp.set_u(f_x * xd );

			(*abs_err)[i][0] = Up[offset];
			(*abs_err)[i][1] =  std::fabs(outp.get_u()- inp.get_u());
			
			residual += std::pow((*abs_err)[i][1], 2);
		}
	//std::cout << "count= " << count << std::endl;

	}
	//std::cout << "Up.size= " << Up.size() << std::endl;
	residual = sqrt(residual / (double)nblig_crop);
	return residual;
}
//CP2KB
double prCameraModelConvert::KBreverseConvert(prPolyCart *input, prCameraModel *output, double theoreticalFOV, vpMatrix *abs_err) 
{
	double residual = 0;
	double nblig = 1024;//(162-53)*2; // image size
	vpColVector x(nblig);
	vpColVector Up(nblig);

	double up = - nblig/2; //u-u0
	double a_0 = input->geta0();
	double a_2 = input->geta2();
	double a_3 = input->geta3();
	double a_4 = input->geta4();

	for (unsigned int i = 0; i < nblig; i++, up++) {
		if (up == 0) {up++;}
		double rho_up_2 = std::pow(std::abs(up), 2);
		double rho_up_3 = std::pow(std::abs(up), 3);
		double rho_up_4 = std::pow(std::abs(up), 4);
		double X_2 = up*up ;
		double Z = (a_0 + a_2 * rho_up_2 + a_3 * rho_up_3 + a_4 * rho_up_4);
		double Z_2 = Z*Z ;
		x[i] =  Z / std::sqrt(X_2+Z_2);
		Up[i] = up;
	}

	vpMatrix A(nblig, 5);
	vpColVector b(nblig);
	vpColVector theta(nblig);

	for (unsigned int i = 0; i < nblig; i++) {

		double theta_i 	= std::acos(x[i]);
		double theta3_i = std::pow(theta_i, 3);
		double theta5_i = std::pow(theta_i, 5);
		double theta7_i = std::pow(theta_i, 7);
		double theta9_i = std::pow(theta_i, 9);

		A[i][0] = theta_i;
		A[i][1] = theta3_i;
		A[i][2] = theta5_i;
		A[i][3] = theta7_i;
		A[i][4] = theta9_i;

		b[i] = ( Up[i] );

		theta[i]=theta_i;
	}
	
	vpColVector params =  A.pseudoInverse() * b;
	//on fixe k_K1 = 1
	double fx = 100;
	double k_k1 = params[0]/fx;
	double k_k2 = params[1]/fx;
	double k_k3 = params[2]/fx;
	double k_k4 = params[3]/fx;
	double k_k5 = params[4]/fx;
	
	std::cout << "fx   = " << fx << "\n";
	std::cout << "k_k1 =" << k_k1 << "\n";
	std::cout << "k_k2 = " << k_k2 << "\n";
	std::cout << "k_k3 = " << k_k3 << "\n";
	std::cout << "k_k4 = " << k_k4 << "\n";
	std::cout << "k_k5 = " << k_k5 << "\n";

	if (abs_err == NULL) {
		for (unsigned int i = 0; i < nblig; i++) {
			prPointFeature inp, outp;
			inp.set_u(Up[i]);
			double Xs,Zs ; 
			double up = Up[i];
			double rho_up_2 = std::pow(std::abs(up), 2);
			double rho_up_3 = std::pow(std::abs(up), 3);
			double rho_up_4 = std::pow(std::abs(up), 4);
			Xs = Up[i];
			Zs = a_0 + a_2 * rho_up_2 + a_3 * rho_up_3 + a_4 * rho_up_4;
	
			double r, xd;
			double theta_i 	= std::acos(Zs/std::sqrt(Zs*Zs+Xs*Xs));
			double theta3_i = std::pow(theta_i, 3);
			double theta5_i = std::pow(theta_i, 5);
			double theta7_i = std::pow(theta_i, 7);
			double theta9_i = std::pow(theta_i, 9);

			r = k_k1*theta_i + k_k2*theta3_i + k_k3*theta5_i + k_k4*theta7_i + k_k5*theta9_i;

			xd = r ;

			outp.set_u(fx * xd );

			residual += std::pow(outp.get_u() - inp.get_u(), 2);
		}
	} else {
		abs_err->resize(nblig, 2, false);
		for (unsigned int i = 0; i < nblig; i++) {
			prPointFeature inp, outp;
			inp.set_u(Up[i]);
			double Xs,Zs ; 
			double up = Up[i];
			double rho_up_2 = std::pow(std::abs(up), 2);
			double rho_up_3 = std::pow(std::abs(up), 3);
			double rho_up_4 = std::pow(std::abs(up), 4);
			Xs = Up[i];
			Zs = a_0 + a_2 * rho_up_2 + a_3 * rho_up_3 + a_4 * rho_up_4;
	
			//double x_FCV = Xs/Zs;
			double r, xd;
			double theta_i 	= std::acos(Zs/std::sqrt(Zs*Zs+Xs*Xs));
			double theta3_i = std::pow(theta_i, 3);
			double theta5_i = std::pow(theta_i, 5);
			double theta7_i = std::pow(theta_i, 7);
			double theta9_i = std::pow(theta_i, 9);

			r = k_k1*theta_i + k_k2*theta3_i + k_k3*theta5_i + k_k4*theta7_i + k_k5*theta9_i;

			xd = r ;

			outp.set_u(fx * xd );

			(*abs_err)[i][0] = Up[i];
			(*abs_err)[i][1] = outp.get_u();//std::fabs(outp.get_u()- inp.get_u());

			residual += std::pow((*abs_err)[i][1], 2);
		}
	}
	residual = sqrt(residual / (double)nblig);
	return residual;
}

//SF2UC
double prCameraModelConvert::convert(prFisheyeEquisolid &fescam, prOmni &ocam, double theoreticalFOV, vpMatrix *abs_err)
{
	double residual = 0;
	unsigned int nblig = vpMath::round(theoreticalFOV)+1; //+1 for zero deg
	vpMatrix A(nblig, 2);
	vpColVector b(nblig);
	double incr_phi = M_PI / 180.; // 1 degree to radian
	double phi = -theoreticalFOV*0.5 * incr_phi;
	
	double fes_au = fescam.getau();
	
	vpColVector Zs(nblig);
	vpColVector Xs(nblig);
	vpColVector ufes(nblig);

	for(unsigned int i = 0 ; i < nblig; i++, phi += incr_phi)
	{
		Zs[i] = cos(phi);
		Xs[i] = sin(phi);
		ufes[i] = 2*fes_au*sin(phi*0.5);
	
		A[i][0] = Xs[i] / ufes[i]; A[i][1] = -1.;
		b[i] = Zs[i];
	}
	
	vpColVector params = A.pseudoInverse()*b;
	
	ocam = prOmni(params[0], params[0], fescam.getu0(), fescam.getv0(), params[1]);

	if(abs_err == NULL)
	{
		//compute residual: later use prCmp?
		for(unsigned int i = 0 ; i < nblig; i++)
		{
			//TODO: use ocam.project3DImage(...)
			residual += pow((params[0]*Xs[i]/(Zs[i]+params[1]) - ufes[i]), 2);
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
			(*abs_err)[i][1] = fabs(params[0]*Xs[i]/(Zs[i]+params[1]) - ufes[i]);
			residual += pow((*abs_err)[i][1], 2);
		}
	}
	
	residual =sqrt(residual)/(double)nblig;
	
	return residual;
}
