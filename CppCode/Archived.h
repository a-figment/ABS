//==== NOT TO BE INCLUDED ====
//@brief Old functions which have not been used for a while
//


// Used for SM^2 conference 
template <typename T>
void WriteTrajectory_Sample(Billiard<T> &Table, const Simulation<T> &Sim)
{

    uint one_gb_in_bytes = 1000000000;
    if (sizeof(T) * 3 * Sim.Time.size() * Sim.nInitialConditions > one_gb_in_bytes)
        std::cout << "WARNING: Trajectory write estimate is greater than 1GB. " << std::endl;

    //int M = Sim.nInitialConditions;
    //if (M > 10)
    //    M = 10;

    // size_t pps = 100;

    // // S
    // std::vector<T> Sx = {0, 0.2, 0.2, 0, 0, 0.2};
    // std::vector<T> Sy = {0, 0, 0.1, 0.1, 0.2, 0.2};
    // std::vector<T> SX;
    // std::vector<T> SY;
    // for (size_t i = 0; i < Sx.size() - 1; ++i)
    // {
    //     std::vector<T> xseq = seq<T>(Sx[i], Sx[i + 1], pps);
    //     std::vector<T> yseq = seq<T>(Sy[i], Sy[i + 1], pps);
    //     for (size_t j = 0; j < xseq.size(); ++j)
    //     {
    //         SX.push_back(xseq[j]);
    //         SY.push_back(yseq[j]);
    //     }
    // }

    // // M
    // std::vector<T> Mx = {0.25, 0.25, 0.35, 0.45, 0.45};
    // std::vector<T> My = {0, 0.2, 0, 0.2, 0};
    // std::vector<T> MX;
    // std::vector<T> MY;
    // for (size_t i = 0; i < Mx.size() - 1; ++i)
    // {
    //     std::vector<T> xseq = seq<T>(Mx[i], Mx[i + 1], pps);
    //     std::vector<T> yseq = seq<T>(My[i], My[i + 1], pps);
    //     for (size_t j = 0; j < xseq.size(); ++j)
    //     {
    //         MX.push_back(xseq[j]);
    //         MY.push_back(yseq[j]);
    //     }
    // }

    // // 2
    // std::vector<T> Tx = {0.475, 0.495, 0.495, 0.475, 0.475, 0.495};
    // std::vector<T> Ty = {0.225, 0.225, 0.2, 0.2, 0.175, 0.175};
    // std::vector<T> TX;
    // std::vector<T> TY;
    // for (size_t i = 0; i < Tx.size() - 1; ++i)
    // {
    //     std::vector<T> xseq = seq<T>(Tx[i], Tx[i + 1], pps);
    //     std::vector<T> yseq = seq<T>(Ty[i], Ty[i + 1], pps);
    //     for (size_t j = 0; j < xseq.size(); ++j)
    //     {
    //         TX.push_back(xseq[j]);
    //         TY.push_back(yseq[j]);
    //     }
    // }

    T d = Table.Obstacles[0]->getVertices()[5].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T dy = Table.Obstacles[0]->getVertices()[1].y;

    std::cout << "dy + d/2 = " << dy + d/2 << std::endl;

    // Point<T> Vel(cos(boost::math::constants::pi<T>() / (4 / (T)1)), sin(boost::math::constants::pi<T>() / (4 / (T)1)));
    // for (int i = 0; i < (SX.size() + MX.size() + TX.size()); ++i)
    //
    // For parallel channel
    //std::vector<T> S = seq<T>(0, 1, Sim.nInitialConditions);
    //Beam<T> B(Point<T>(-0.45, d / 4), Point<T>(-0.2, dy + 0.1), 0, 0);
    
    // For Lorentz gas
    std::vector<T> S = seq<T>(0, 1, Sim.nInitialConditions);
    Beam<T> B(Point<T>(0.01, 0.25), Point<T>(0.25, 0.01), 0, 0);

    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        // Particle<T> P;

        // if (i < SX.size())
        //     P = Particle<T>(Point<T>(SX[i], SY[i]), Vel / Norm<T>(Vel), -1, -1, Sim.Time.back());
        // if (i >= SX.size() && i < (MX.size() + SX.size()))
        //     P = Particle<T>(Point<T>(MX[i - SX.size()], MY[i - SY.size()]), Vel / Norm<T>(Vel), -1, -1, Sim.Time.back());
        // if (i >= SX.size() + MX.size() && i < (MX.size() + SX.size() + TX.size()))
        //     P = Particle<T>(Point<T>(TX[i - SX.size() - MX.size()], TY[i - SY.size() - MY.size()]), Vel / Norm<T>(Vel), -1, -1, Sim.Time.back());
        // P.Diameter = Point<T>(1, 1);

        Particle<T> P = Table.getInitialState(0, -1, Sim.Time.back());
	P.Pos = B.getPointOnBeam(S[i]);
        //P.Vel = Point<T>(cos(-25 * 3.14/180), sin(-25 * 3.14 / 180));
        P.Vel = Point<T>(cos(43 * 3.14/180), sin(43 * 3.14 / 180));
	
        //Particle<T> PC = P;
	//
        //Particle<T> P = Table.getIC_FixedLine(0, -1, Sim.TotalTime, Point<T>(-0.45, d / 4), Point<T>(-0.2, dy + 0.1));

        //Particle<T> P = Table.getIC_ExpandingBall(0, -1, Sim.TotalTime, Point<T>(dx/2, (-dy/dx)*((dx/2) - dx) + d/2), 0.001);

        //Trajectory<T> Traj = Table.getTrajectorySample(P, Sim.Time); // Continuous and distrete time (for plotting)
        //Traj.Write(Sim.DataDestination + "Traj_" + std::to_string(i) + ".dat");

        Trajectory<T> ContTraj = Table.getTrajectory(P, Sim.Time); // Continuous time only
        if (Table.Failed)
            std::cout << "Failed" << std::endl;
        ContTraj.WritePositions(Sim.DataDestination + "Positions_" + std::to_string(i) + ".dat");

        //Trajectory<T> DiscTraj = Table.getTrajectory(PC, Sim.nCollisions); // Collisions only
        //DiscTraj.Write(Sim.DataDestination + "DiscTraj_" + std::to_string(i) + ".dat");
    }
}


// // @brief Trajectories are conditionally generated, requiring that |x(t)| > Sim.L.back()
// // @Write - All the usual metrics of a finite system
// template <typename T>
// inline void EscapeToSymbolicSequence(Billiard<T> &Table, Simulation<T> &Sim)
// {
//     Ensemble<T> MeanEscapeTime(Sim.Symbols.size());
//     std::vector<Ensemble<T>> FiniteMoments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));

//     std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//     for (int i = 0; i < Sim.nInitialConditions; ++i)
//     {
//         if (10 * i % Sim.nInitialConditions == 0 && i != 0)
//             std::cout << "Simulated " << i << " initial conditions." << std::endl;

//         // For finite trajectories we need discrete trajectories
//         EscapeData<T> ED(Sim.Symbols.size());

//         Particle<T> IC = Table.getInitialState(0, -1, Sim.TotalTime);
//         Trajectory<T> DiscreteTraj = Table.getFiniteTrajectory(IC, Sim.nCollisions);

//         if (Table.Failed)
//         {
//             Table.ErrorReport(IC, i);
//             i--; // Want M initial conditions
//         }

//         ED = DiscreteTraj.getEscapeData(Sim.Symbols);
//         MeanEscapeTime.Push(ED.EscapeTimes);

//         // Finite moments
//         Trajectory<T> FiniteContinuousTraj = DiscreteTraj;
//         FiniteContinuousTraj.MakeContinuous(Sim.Time);
//         std::vector<std::vector<T>> DXp = FiniteContinuousTraj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
//         for (int j = 0; j < Sim.nMoments; ++j)
//             FiniteMoments[j].Push(DXp[j]);

//         WriteVector(Sim.DataDestination + "FP_Time.dat", ED.EscapeTimes);
//         WriteVector(Sim.DataDestination + "FP_Abscissa.dat", ED.EscapePosition_x);
//         WriteVector(Sim.DataDestination + "FP_Ordinate.dat", ED.EscapePosition_y);
//         WriteVector(Sim.DataDestination + "FP_Direction.dat", ED.EscapeDirection);
//         DiscreteTraj.getParticle(0).WriteState(Sim.DataDestination + "InitialConditions.dat");
//     }
//     WriteVector(Sim.DataDestination + "MFPT.dat", MeanEscapeTime.PullMeans());
//     for (size_t i = 0; i < FiniteMoments.size(); ++i)
//         FiniteMoments[i].Write(Sim.DataDestination + "FiniteMoment_" + std::to_string(i) + ".dat");

//     std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//     const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
//     Sim.Print(RunTime);
//     Sim.Write();
// }


//template <typename T>
//inline void getDisplacementSummary_t(Billiard<T> &Table, Simulation<T> &Sim)
//{
//    // Stream for writing displacements
//    std::ofstream Writer_DX(Sim.DataDestination + "DX.dat", std::ios::out | std::ios::app);
//    Writer_DX << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
//
//    // Moments of displacement
//    std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));
//
//    // Displacement correlatoin
//    //std::vector<size_t> h = { 0, 10, 50, 100, 500, 1000, 5000 }; // Assumes time step of unity
//    //std::vector<size_t> h = { 0, 50, 100, 500, 1000, 5000 }; // Assumes time step of unity
//    //std::vector<size_t> h = { 0, 1, 2, 10, 20 }; // Assumes time step of 1000
//    std::vector<size_t> h = { 1, 2, 10, 20, 100 }; // Assumes time step of 10
//    std::vector<Ensemble<T>> PAC(h.size(), Ensemble<T>(Sim.Time.size() - h[0]));
//
//    std::vector<Ensemble<T>> TAMSD(h.size(), Ensemble<T>(Sim.Time.size() - h[0]));
//
//    std::vector<Ensemble<T>> PACL(h.size(), Ensemble<T>(Sim.Time.size() - h[0]));
//    std::vector<Ensemble<T>> PACR(h.size(), Ensemble<T>(Sim.Time.size() - h[0]));
//
//    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
//    std::cout << "Indices found: ";
//    for (size_t i = 0; i < Indices.size(); ++i)
//        std::cout << Indices[i] << " ";
//    std::cout << std::endl;
//    std::cout << "Collecting displacements at times t = ";
//    for (size_t i = 0; i < Indices.size(); ++i)
//        std::cout << Sim.Time[Indices[i]] << " ";
//    std::cout << std::endl;
//
//    // Main loop
//    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//    for (int i = 0; i < Sim.nInitialConditions; ++i)
//    {
//        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
//            std::cout << "Simulated " << i << " initial conditions." << std::endl;
//
//        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);
//
//        // Extract trajectory properties
//        std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
//        for (int j = 0; j < Sim.nMoments; ++j)
//            Moments[j].Push(DXp[j]);
//
//        // Get PAC
//        std::vector<T> xDis = Traj.getDisplacements_x();
//        for (size_t i = 0; i < h.size(); ++i) {
//            std::vector<T> pac_temp(xDis.size() - h[i]); // trailing zeros
//
//            std::vector<T> pac_temp_L(xDis.size() - h[i]); // trailing zeros
//            std::vector<T> pac_temp_R(xDis.size() - h[i]); // trailing zeros
//            for (size_t j = 0; j < xDis.size() - h[i]; ++j) {
//                pac_temp[j] = xDis[j] * xDis[j + h[i]];
//
//		if (((xDis[j] / Sim.Time[j]) < 0.4) && ((xDis[j + h[i]] / Sim.Time[j + h[i]]) < 0.4)) {
//		    pac_temp_L[j] = xDis[j] * xDis[j + h[i]];
//		} else {
//		    pac_temp_R[j] = xDis[j] * xDis[j + h[i]];
//		}
//            }
//            PAC[i].Push(pac_temp);
//            PACL[i].Push(pac_temp_L);
//            PACR[i].Push(pac_temp_R);
//        }
//
//        // Set and write displacements
//        for (size_t j = 0; j < Indices.size(); ++j) 
//            Writer_DX << xDis[Indices[j]] << " ";
//        Writer_DX << std::endl;
//    }
//
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
//
//    for (size_t i = 0; i < Moments.size(); ++i)
//        Moments[i].Write(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat");
//
//    for (size_t i = 0; i < PAC.size(); ++i) {
//        PAC[i].Write(Sim.DataDestination + "PAC_" + std::to_string(i) + ".dat");
//        PACL[i].Write(Sim.DataDestination + "PACL_" + std::to_string(i) + ".dat");
//        PACR[i].Write(Sim.DataDestination + "PACR_" + std::to_string(i) + ".dat");
//    }
//    Sim.Print(RunTime);
//}


//template <typename T>
//inline void ContinuousSystemSummary(Billiard<T> &Table, Simulation<T> &Sim)
//{
//    std::ofstream Writer_DX(Sim.DataDestination + "DX.dat", std::ios::out | std::ios::app);
//    Writer_DX << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
//
//    std::vector<size_t> h = { 0, 10, 50, 100, 500, 1000, 5000 }; // Assumes time step of unity
//
//    // Involves all initial conditions
//    //std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));
//    std::vector<Ensemble<T>> PAC(h.size(), Ensemble<T>(Sim.Time.size() - h[1]));
//
//    //std::vector<int> Indices = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
//    std::vector<int> Indices = { 601, 881, 942, 979, 1004 };
//
//    std::cout << "Collecting displacements at times t = ";
//    for (size_t i = 0; i < Indices.size(); ++i)
//        std::cout << Sim.Time[Indices[i] - 1] << " ";
//    std::cout << std::endl;
//
//    std::vector<T> xDisplacement(Indices.size(), -1);
//
//    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//    for (int i = 0; i < Sim.nInitialConditions; ++i)
//    {
//        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
//            std::cout << "Simulated " << i << " initial conditions." << std::endl;
//
//        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);
//
//        // Extract trajectory properties
//        //std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
//        //for (int j = 0; j < Sim.nMoments; ++j)
//        //    Moments[j].Push(DXp[j]);
//
//        // Get 
//        std::vector<T> xDis = Traj.getDisplacements_x();
//        for (size_t i = 0; i < h.size(); ++i) {
//            std::vector<T> pac_temp(xDis.size() - h[i]); // trailing zeros
//            for (size_t j = 0; j < xDis.size() - h[i]; ++j) {
//                pac_temp[j] = xDis[j] * xDis[j + h[i]];
//            }
//            PAC[i].Push(pac_temp);
//        }
//
//        // Set and write displacements
//        for (size_t j = 0; j < xDisplacement.size(); ++j)
//            xDisplacement[j] = xDis[Indices[j] - 1];
//        for (size_t i = 0; i < xDisplacement.size(); ++i)
//            Writer_DX << xDisplacement[i] << " ";
//        Writer_DX << std::endl;
//
//        //WriteVector<T>(Sim.DataDestination + "xPositions.dat", xPositions);
//        //WriteVector<T>(Sim.DataDestination + "xDisplacements.dat", xDisplacement);
//        //Traj.MakeParticle(0).WriteState(Sim.DataDestination + "InitialConditions.dat");
//    }
//
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
//
//    //for (size_t i = 0; i < Moments.size(); ++i)
//    //    Moments[i].Write(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat");
//    for (size_t i = 0; i < PAC.size(); ++i)
//        PAC[i].Write(Sim.DataDestination + "PAC_" + std::to_string(i) + ".dat");
//    Sim.Print(RunTime);
//}


// @brief Writes displacements at times t[Indices] 
template <typename T>
inline void getFront(Billiard<T> &Table, Simulation<T> &Sim, const T& vbal)
{
    std::ofstream Writer_DX(Sim.DataDestination + "DX.dat", std::ios::out | std::ios::app);
    Writer_DX << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    std::ofstream Writer_FP(Sim.DataDestination + "FP.dat", std::ios::out | std::ios::app);
    std::ofstream Writer_IC(Sim.DataDestination + "IC.dat", std::ios::out | std::ios::app);
    //std::ofstream Writer_TSD(Sim.DataDestination + "TSD.dat", std::ios::out | std::ios::app);

    // Geometry 
    //T d  = Table.Obstacles[0]->getVertices()[5].y;
    //T dy = Table.Obstacles[0]->getVertices()[1].y;
    //T dx = Table.Obstacles[0]->getVertices()[0].x;

    T eps = 0.01;
    //T eps = 0.055;
    T Vbal = vbal;

    // Some angles
    //std::vector<T> th = seq<T>(0.2, 2.5 / 4, Sim.nInitialConditions);
    //std::vector<size_t> nCollisions(Sim.nInitialConditions, 0);
    
   
    // Get the mean value of dx / t for each particle (time average)
    OnlineStatistics<T> OS; // For finding the average
    std::vector<T> VB_TimeAvg(Sim.nInitialConditions);

    // Get the ensemble average of dx/t for each t (ensemble average)
    Ensemble<T> VB_EnseAvg(Sim.Time.size());


    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = ";
    for (size_t i = 0; i < Indices.size(); ++i)
        std::cout << Sim.Time[Indices[i]] << ", ";
    std::cout << " with indices = ";
    for (size_t i = 0; i < Indices.size(); ++i)
        std::cout << Indices[i] << ", ";
    std::cout << std::endl;

    size_t M = Sim.nInitialConditions; //
    size_t count = 0;
    bool flag = false;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < M; ++i)
    {
	bool ThreshCond = true;
        if (10 * count % Sim.nInitialConditions == 0 && count != 0 && flag) {
            std::cout << "Simulated " << count << " initial conditions." << std::endl;
	    flag = false;
	}

        // ON P.B.
        //Particle<T> P0 = Table.getIC_FixedLine(0, -1, Sim.TotalTime, Point<T>(-dx, (T)0), Point<T>(-dx, d));
        //P0.Vel = Point<T>(1,0);
        //
        //Particle<T> P0 = Table.getIC_FixedPosition(0, -1, Sim.TotalTime, Point<T>(0, dy + 0.0001));
        //P0.Vel = Point<T>(cos(th[i]), sin(th[i]));
        //std::cout << P0 << std::endl;
        
        // Standard distribution
        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Particle<T> Pn = P0;

        Trajectory<T> Traj = Table.getTrajectory(Pn, Sim.Time);

        std::vector<T> xDis = Traj.getDisplacements_x();

        // Rejection sample
        for (size_t j = 0; j < Indices.size(); ++j) {
            T v = abs(xDis[Indices[j]]) / Sim.Time[Indices[j]]; // Make sure v = .5 at all times

	    if (v < Vbal - eps || v > Vbal + eps) { // Only one instance of dx/t < .2 removes the sample
		ThreshCond = false;
		break;
	    }
	}

	if (ThreshCond) {
	    //std::cout << "FOUND!" << " M = " << M << ", i = " << i << std::endl;
	    for (size_t j = 0; j < Indices.size(); ++j) 
		Writer_DX << xDis[Indices[j]] << " ";
	    Writer_DX << std::endl;


	    // Rescale dx / t, take time and ensemble average
	    //std::vector<T> TSD_DX(xDis.size());
	    //for (size_t k = 0; k < xDis.size(); ++k) {
	    //   TSD_DX[k] = xDis[k] / Sim.Time[k];
	    //   OS.Update(TSD_DX[k]);
	    //}
	    //VB_TimeAvg[i] = OS.getMean();
	    //OS.Clear();
	    //VB_EnseAvg.Push(TSD_DX);	

	    // Write dx / t of individual trajectory
	    //WriteVector<T>(TSD_DX, Writer_TSD);

	    // Write initial conditions
            WriteVector<T>({P0.Pos.x, P0.Pos.y, P0.Vel.x, P0.Vel.y}, Writer_IC);
	    flag = true;
	    count++;


	    // Need to get collisions until some maximum time is reached and then determing

	    // 
	    Particle<T> Pc = P0, Pt = P0;
            Trajectory<T> Traj1 = Table.getTrajectory(Pc, Sim.nCollisions);
            //Trajectory<T> Traj1 = Table.getTrajectory(Pc, Sim.nCollisions, Sim.Time[Indices[0]], Sim.Time[Indices.back()-1]);
	    
	    std::vector<T> LFP = Traj1.getLiftedFreePaths();
            WriteVector<T>(LFP, Writer_FP);

	    ///Trajectory<T> Traj2 = Table.getTrajectory(Pt, Sim.Time);
            //Traj1.WritePositions(Sim.DataDestination + "ColTraj-" + std::to_string(i) + ".dat");
            //Traj2.WritePositions(Sim.DataDestination + "ConTraj-" + std::to_string(i) + ".dat");

        } else {
	    M++;
	}
    }

    //VB_EnseAvg.Write(Sim.DataDestination + "eMeanDXT.dat");
    //WriteVector(Sim.DataDestination + "tMeanDXT.dat", VB_TimeAvg);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    Sim.Print(RunTime);
}

template <typename T>
void LaminarICSample(Billiard<T>& Table, const std::vector<T>& TIME, 
	const uint& N, 
	const int& M, 
	const T& minimum_step_distance, 
	const std::vector<size_t>& nBB, const std::string& DataDir) {

    // ------------- Fiddling around with IC's ---------------
    T d = Table.Obstacles[0]->getVertices()[5].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T dy = Table.Obstacles[0]->getVertices()[1].y;

    // -------------------------------------------------------

    // --- WRITERS ---
    Writer<size_t> w_DC(DataDir + "DC.dat", true); // Number of direction changes
    std::vector<Writer<T>> w_JUMP(nBB.size());
    for (size_t i = 0; i < nBB.size(); ++i) {
	w_JUMP[i] = Writer<T>(DataDir + "n" + std::to_string(nBB[i]) + "_jump.dat", true);
    }

    T eps = 0.0125;
    //T vb_lb = (1/T(2)) - eps;
    //T vb_ub = (1/T(2)) + eps;
    T vb_lb = 2825/(T)10000 - eps;
    T vb_ub = 2825/(T)10000 + eps;

    std::vector<size_t> nDC(nBB.size(), 0);
    // Main loop
    for (int k = 0; k < M; ++k) {
	if (10 * k % M == 0 && k != 0)
	    std::cout << "Simulated " << k << " initial conditions." << std::endl;

	Particle<T> P = Table.getInitialState(0, -1, TIME.back());
	T theta = atan(P.Vel.y / P.Vel.x);

	// Get trajectory - discrete
	Trajectory<T> Traj = Table.getTrajectory(P, N);
	std::vector<T> x = Traj.getLiftedPositions_x();
	std::vector<Point<T>> Vel = Traj.getVelocities();
	std::vector<T> t = Traj.getLiftedFreePaths();
	std::vector<Point<long int>> cid = Traj.getCellIDs();

	// curloop
	std::vector<int> FlowDirections = getFlowDirections(x, Vel, cid, dy, dx);
	
	for (size_t i = 0; i < nBB.size(); ++i) {

	    std::vector<Point<T>> lambda = getLaminarSteps(t, x, FlowDirections, minimum_step_distance, nBB[i]);
	    //std::vector<Point<T>> lambda = getLaminarStepsFD(t, x, FlowDirections, minimum_step_distance, nBB[i]);
	    nDC[i] += lambda.size();

	    if (lambda.size() > 0) {
		bool found = false;
		for (size_t j = 0; j < lambda.size(); ++j) {
		    T v = std::abs(lambda[j].y) / lambda[j].x;  
		    if (v > vb_lb && v < vb_ub) {
			found = true;
			break;
		    }
		}
		
		if (found) {
		    std::vector<T> ths(lambda.size(),0);
		    WriteColVector_app<T>(DataDir + "theta.dat",  { theta });
		    //for (size_t j = 0; j < lambda.size(); ++j) 
		    //    ths[j] = theta;
		    std::vector<T> lx = getCoordFromVectorPoint(lambda, "x");
		    std::vector<T> ly = getCoordFromVectorPoint(lambda, "y");
		    //w_JUMP[i].WriteColVectors(lx, ly, ths);
		    w_JUMP[i].WriteColVectors(lx, ly);
		}

	    }
	}
    }
    w_DC.WriteRowVector(nDC);
}

template <typename T>
std::vector<Point<T>> LaminarOnlineFilter(const std::vector<size_t>& CC, const size_t& nBB, const std::vector<Point<T>>& steps, const Point<T>& minstep) {

    assert(CC.size() == steps.size());

    // Steps to return
    std::vector<Point<T>> filtered_steps;
    filtered_steps.reserve(steps.size());

    size_t i = 0;

    while (i < CC.size()-1) { // n steps forward, k < n steps back

	size_t start_idx = i;
	size_t j = i;
	while (j < CC.size()-1 && CC[j+1] < CC[j]) 
	    j++;
	size_t end_idx = j;

	Point<T> step_sum = steps[start_idx];
	for (size_t j = start_idx; j < end_idx; ++j)
	    step_sum = step_sum + steps[j];

        if (step_sum.x > minstep.x && std::abs(step_sum.y) > minstep.y)
	   filtered_steps.push_back(step_sum);

	i = end_idx + 1;
    }


    //size_t start_idx = 0;
    //for (size_t i = 0; i < CC.size()-1; i += 1) {

    //    if (CC[i] < CC[i+1]) {

    //        Point<T> step_sum = steps[start_idx];
    //        for (size_t j = start_idx+1; j < i+1; ++j)
    //    	step_sum = step_sum + steps[j];

//        if (step_sum.x > minstep.x && std::abs(step_sum.y) > minstep.y)
//	   filtered_steps.push_back(step_sum);

    //        start_idx = i+1;
    //        i++;
    //    }
    //}


    //size_t start_idx = 0;
    //for (size_t i = 0; i < CC.size()-1; ++i) {

    //    if (CC[i+1] > nBB) {

    //        Point<T> step_sum = steps[start_idx];
    //        for (size_t j = start_idx+1; j < i+1; ++j)
    //    	step_sum = step_sum + steps[j];

//        if (step_sum.x > minstep.x && std::abs(step_sum.y) > minstep.y)
//	   filtered_steps.push_back(step_sum);

    //        start_idx = i+1;
    //    }
    //}

    //for (size_t i = 0; i < steps.size(); ++i) {
    //    T v = std::abs(steps[i].y) / steps[i].x;
    //    if (v < 4/(T)10)
    //        filtered_steps.push_back(steps[i]);
    //}

    //if (steps.size() > 2) {

    //    size_t i = 0;
    //    while (i < CC.size()-2) {

    //        int cumsum = 0;
    //        int start_idx = i;

    //        if (cumsum >= 0) { // Necessary to avoid b2b while loops
    //    	while (i < CC.size()-2 && cumsum >= 0) {
    //    	    cumsum += (CC[i] - CC[i+1]);
    //    	    i += 2;
    //    	}
    //        } else {
    //    	while (i < CC.size()-2 && cumsum < 0) {
    //    	    cumsum += (CC[i] - CC[i+1]);
    //    	    i += 2;
    //    	}
    //        }

    //        Point<T> step_sum = steps[start_idx];
    //        int end_idx = i - 2; // The last update triggered loop condition, backpeddle Point<T> step_sum = steps[start_idx];
    //        for (int j = start_idx+1; j < end_idx; ++j)
    //    	step_sum = step_sum + steps[j];

    //        if (std::abs(step_sum.y) > minstep.y)
    //    	filtered_steps.push_back(step_sum);

    //        i = i + 2; // Resume
    //    }
    //}


    //size_t K = nBB;
    //if (steps.size() >= K) {
    //    for (size_t i = 0; i < CC.size()-K; i += K) {
    //        size_t outer_sum = CC[i] + CC[i + K - 1];
    //        size_t inner_sum = 0;
    //        for (size_t k = i+1; k < i + K - 1; ++k)
    //    	inner_sum += CC[k];
    //        if (outer_sum > inner_sum) { // Condition on if jump should extend 
    //    	Point<T> ss = steps[i+1];
    //    	for (size_t k = i+2; k < i + K - 1; ++k) // SUM of all steps BETWEEN and INCLUDING large
    //    	    ss = ss + steps[k];

    //    	//Point<T> ss = steps[i+1];
    //    	//for (size_t k = i+2; k < i + K - 2; ++k) // The SUM of steps ONLY BETWEEN large steps
    //    	//    ss = ss + steps[k];

    //    	//for (size_t k = i+1; k < i + K - 1; ++k) // The INDIVIDUAL steps between large steps
    //        if (std::abs(step_sum.y) > minimum_step_distance)
    //    	filtered_steps.push_back(step_sum);

    //    	if (std::abs(ss.y) > minimum_step_distance) 
    //    	    filtered_steps.push_back(ss);
    //        }
    //    }
    //}
    return filtered_steps;
}


template <typename T>
void LaminarOnline(Billiard<T>& Table, const std::vector<T>& TIME, 
	const uint& N, 
	const int& M, 
	const Point<T>& minstep, 
	const std::vector<size_t>& nBB, const bool& filter, const size_t& max_num_steps, const std::string& DataDir) {

    // ------------- Fiddling around with IC's ---------------
    T d = Table.Obstacles[0]->getVertices()[5].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T dy = Table.Obstacles[0]->getVertices()[1].y;
    //std::vector<Particle<T>> Ps = generateICs(Table, TIME, N, M, minstep, DataDir);
    //WriteParticles<T>(DataDir + "P0.dat", Ps);

    T mid = PI/300;
    std::vector<T> theta_left = seq<T>(-mid, -1e-5, M/2);
    std::vector<T> theta_right = seq<T>(1e-5, mid, M/2);
    std::vector<T> theta_ics = theta_left;
    theta_ics.insert(theta_ics.end(), theta_right.begin(), theta_right.end());

    // -------------------------------------------------------

    // ------------------ WRITERS ------------------
    size_t nWriters = 1; // Filter is off
    if (filter)
	nWriters = nBB.size();
    std::vector<Writer<T>> w_JUMP(nWriters);
    for (size_t i = 0; i < nWriters; ++i) {
	if (filter)
	    w_JUMP[i] = Writer<T>(DataDir + "N" + std::to_string(nBB[i]) + "_jump.dat", true);
	else 
	    w_JUMP[i] = Writer<T>(DataDir + "N" + std::to_string(0) + "_jump.dat", true);
    }
    //Writer<size_t> w_DC(DataDir + "DC.dat", true); // Number of direction changes
    //Writer<size_t> w_CC(DataDir + "CC.dat", true); // Number of direction changes
    //Writer<size_t> w_TL(DataDir + "TL.dat", true); // Trajectory length
    // ---------------------------------------------

    std::vector<size_t> ncoll = { 1000, 10000, 100000, 1000000 };
    std::vector<std::vector<T>> TrajLength(ncoll.size(), std::vector<T>(M,0));
    std::vector<size_t> nDC(nWriters, 0); // Number of rows in file

    // Main loop
    for (int k = 0; k < M; ++k) {
	if (10 * k % M == 0 && k != 0)
	    std::cout << "Simulated " << k << " initial conditions." << std::endl;

	//Particle<T> P = Ps[k];
	Particle<T> P = Table.getInitialState(0, -1, TIME.back());
	//P.Pos = Point<T>(0, (d + dy)/2);
	//P.Pos = Point<T>(0, 3*d/(T)2);
	//P.Vel = Point<T>(cos(theta_ics[k]), sin(theta_ics[k]));
	//T th = theta_ics[k];

	Point<T> diam = Table.getDiameter();
	std::vector<Point<T>> steps;
	std::vector<size_t> CC;

	T x_init = P.LPos(diam).x;
	//steps.push_back(Point<T>(0, x_init)); // Needed for reconstruction
	int start_flow = getFlowDirection(x_init, P.Vel, P.CellID, dy, dx); // initial flow direction

	//std::vector<T> x, t, y, vx, vy;
	//x.push_back(x_init);
	//t.push_back(0);
	//y.push_back(P.Pos.y);

	int nDirChange = 0;
	int dir = start_flow;
	size_t n = 0, tl_idx = 0;
	while (n < N) {
	//while (P.TotalTime < TIME.back()) {
	    T lt = 0, lx = 0, lx_init = x_init;
	    int steps_in_flow = 0; // Num steps along direction of flow
				   
	    while (start_flow == dir) {
		x_init = P.LPos(diam).x; // Position 
		T tt = P.TotalTime;
		Table.Update(P);
		//x.push_back(P.LPos(diam).x);
		//t.push_back(P.TotalTime);
		//y.push_back(P.Pos.y);

		lt += (P.TotalTime - tt); 
		start_flow = getFlowDirection(P.LPos(diam).x, P.Vel, P.CellID, dy, dx);
		steps_in_flow++;
		n++;

		if (n == ncoll[tl_idx]) {
		    TrajLength[tl_idx][k] = P.TotalTime;
		    tl_idx++;
		}

	    }
	    lx = P.LPos(diam).x - lx_init;
	    x_init = P.LPos(diam).x;

	    if (!filter) { // Jump events to be written
			   //if (P.TotalTime < TIME.back()) {
		if (lt > minstep.x && std::abs(lx) > minstep.y) { 
		    steps.push_back(Point<T>(lt, lx));
		    CC.push_back(steps_in_flow); 
		} 
		//}
	    } else { // When filtering, store all jump events
		steps.push_back(Point<T>(lt, lx));
		CC.push_back(steps_in_flow); 
	    }
	    dir = -dir;

	    nDirChange++;
	}

	//std::cout << P.TotalTime << " " << steps.size() << std::endl;
	if (!filter) {
	    if (steps.size() > 0) { 
		nDC[0] += steps.size();
		if (nDC[0] < max_num_steps) {
		    std::vector<T> lx = getCoordFromVectorPoint(steps, "x");
		    std::vector<T> ly = getCoordFromVectorPoint(steps, "y");
		    //std::vector<T> ths(lx.size(),0);
		    //std::vector<int> ndir(lx.size(),0);
		    //for (size_t i = 0; i < lx.size(); ++i) {
		    //    ths[i] = th;
		    //    ndir[i] = nDirChange;
		    //}
		    //w_JUMP[0].WriteColVectors(lx, ly, ths, ndir);
		    w_JUMP[0].WriteColVectors(lx, ly, CC);
		} else {
		    return; // Nothing else to write, so exit
		}
	    }
	} else {
	    for (size_t j = 0; j < nBB.size(); ++j) { 
		std::vector<Point<T>> filtered_steps = LaminarOnlineFilter(CC, nBB[j], steps, minstep);
		nDC[j] += filtered_steps.size();
		std::vector<T> lx = getCoordFromVectorPoint(filtered_steps, "x");
		std::vector<T> ly = getCoordFromVectorPoint(filtered_steps, "y");

		if (filtered_steps.size() > 0 && nDC[j] < max_num_steps)
		    w_JUMP[j].WriteColVectors(lx, ly);//#, tt, th);

	    }
	}
	//WriteColVector<T>(DataDir + "n" + std::to_string(nBB[0]) + "_" + std::to_string(k) + "--t.dat", t);
	//WriteColVector<T>(DataDir + "n" + std::to_string(nBB[0]) + "_" + std::to_string(k) + "--x.dat", x);
	//WriteColVector<T>(DataDir + "n" + std::to_string(nBB[0]) + "_" + std::to_string(k) + "--y.dat", y);
	//w_CC.WriteRowVector(CC);
    }
//    for (size_t i = 0; i < TrajLength.size(); ++i)
//	WriteColVector<T>(DataDir + std::to_string(i) + "_TL.dat", TrajLength[i]);
}



// @brief Consider long trajectories, and compare the distribution of jumps for a filter applied to
//	  smaller subsets of that trajectory. 
template <typename T>
void JumpScaling(Billiard<T>& Table, const std::vector<T>& TIME, 
	const std::vector<size_t>& Ti, 
	const std::vector<uint>& W, 
	const int& M, 
	const T& dxs, 
	const T& dxt, 
	const T& mean_pd, const std::string& DataDir) {

    // Count the total number of samples
    //std::vector<size_t> K(W.size(),0);
    std::vector<std::vector<size_t>> K(W.size(),std::vector<size_t>(Ti.size(), 0));

    // Number of trajectories which failed graining for each window
    std::vector<uint> ungrainable(W.size(), 0);

    // --- WRITERS ---
    //std::vector<std::vector<Writer<T>>> w_JUMP(W.size(), std::vector<Writer<T>>(Ti.size()));
    //for (size_t i = 0; i < W.size(); ++i) {
    //    std::vector<Writer<T>> w_temp(Ti.size());
    //    for (size_t j = 0; j < Ti.size(); ++j) {
    //        std::string file = DataDir + "w" + std::to_string(W[i]) + "_T" + std::to_string(j) + "_jump.dat";
    //        w_JUMP[i][j] = Writer<T>(file, true);
    //        //
    //        //temp[j] = Writer<T>(file, true);
    //    }
    //    //w_JUMP.push_back(w_temp);
    //}


    //std::vector<Writer<T>> w_JUMP(W.size());
    //for (size_t i = 0; i < W.size(); ++i) {
    //    w_JUMP[i] = Writer<T>(DataDir + "w" + std::to_string(W[i]) + "_jump.dat", true);
    //}

    // Main loop
    for (int k = 0; k < M; ++k) {
        if (10 * k % M == 0 && k != 0)
            std::cout << "Simulated " << k << " initial conditions." << std::endl;

        Particle<T> P = Table.getInitialState(0, -1, TIME.back());
	// Get trajectory - continuous
    	Trajectory<T> Traj = Table.getTrajectory(P, TIME);
    	std::vector<T> x = Traj.getLiftedPositions_x();
	//WriteColVector<T>(DataDir + std::to_string(k) + "_x.dat", x);

	// Apply window filter to current particle position x(t)
	for (size_t i = 0; i < W.size(); ++i) {

	    for (size_t j = 0; j < Ti.size(); ++j) { 

		std::vector<T> new_time = index(TIME, 0, Ti[j]);
		std::vector<T> new_pos = index(x, 0, Ti[j]);
		CoarseGrain<T> cg(new_time, new_pos);

		cg.Grain(W[i]);
		if (cg.isGrainable()) { // 

		    JumpEvent<T> JE = cg.getJumpEvent(dxt, dxs, mean_pd);

		    // Writing raw jump time and jump distance
		    K[i][j] += JE.count.size(); 
		    if (K[i][j] < (size_t)1e6 && JE.count.size() > 0) {
			//w_JUMP[i][j].WriteColVector(JE);

			//JE.Write(DataDir + "w" + std::to_string(W[i]) + "_T" + std::to_string(j) + "_jump.dat");
			WriteColVector<T>(DataDir + "w" + std::to_string(W[i]) + "_T" + std::to_string(j) + "_jump.dat", JE);
			//WriteColVector<T>(DataDir + "w" + std::to_string(W[i]) + "_" + std::to_string(k) + "--x.dat", x);
			//WriteColVector<T>(DataDir + "w" + std::to_string(W[i]) + "_" + std::to_string(k) + "--t.dat", t);
		    }

		} else {
		    ungrainable[i]++;
		}
	    }
	} 
    }
}

// Allows for @nStepsBack 
// @return [vector{point}] Time (x) and step distance (y)
template <typename T>
std::vector<Point<T>> getLaminarSteps(const std::vector<T>& t, const std::vector<T>& x, const std::vector<int>& FD, const Point<T>& minstep, const size_t& nStepsBack) {

    std::vector<Point<T>> steps;
    steps.reserve(FD.size());
    std::vector<size_t> CC = ConsecutiveCounts(FD);

    //std::cout << FD << std::endl;
    //std::cout << CC << std::endl;
    //std::cout << t << std::endl;
    //std::cout << x << std::endl;

    steps.push_back(Point<T>(t[0], x[0]));

    std::vector<size_t> ccc(CC.size(), 0);
    std::partial_sum(CC.begin(), CC.end(), ccc.begin());
    ccc.insert(ccc.begin(), 0);

    size_t idx = 0;
    for (size_t i = 0; i < ccc.size()-1; ++i) {

	size_t next_step_length = ccc[i+1] - ccc[i];

	//if (next_step_length > nStepsBack) {
	if (next_step_length > nStepsBack) {
	    T lt = 0;
	    for (size_t j = ccc[idx]+1; j <= ccc[i+1]; ++j) 
	        lt += t[j];
	    T lx = x[ccc[i+1]] - x[ccc[idx]];

	    if (lt > minstep.x && std::abs(lx) > minstep.y)
	        steps.push_back(Point<T>(lt, lx));
	    idx = i+1;
	}
    }
    return steps;
}

// @brief Computes laminar stretches - jump time and displacement
//
template <typename T>
void Laminar(Billiard<T>& Table, const std::vector<T>& TIME, 
	const uint& N, 
	const int& M, 
	const Point<T>& minstep, 
	const std::vector<size_t>& nBB, const std::string& DataDir) {

    // ------------- Fiddling around with IC's ---------------
    //T d = Table.Obstacles[0]->getVertices()[5].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T dy = Table.Obstacles[0]->getVertices()[1].y;

    //T eps = (0.05)*PI/(T)180;
    //std::vector<T> thetas = seq<T>(-eps, eps, M);               
    // -------------------------------------------------------

    // --- WRITERS ---
    //std::vector<Writer<T>> w_VEL(nBB.size());
    std::vector<Writer<T>> w_JUMP(nBB.size());
    for (size_t i = 0; i < nBB.size(); ++i) {
	w_JUMP[i] = Writer<T>(DataDir + "n" + std::to_string(nBB[i]) + "_jump.dat", true);
	//w_VEL[i] = Writer<T>(DataDir + "n" + std::to_string(nBB[i]) + "_vel.dat", true);
    }
    Writer<size_t> w_DC(DataDir + "DC.dat", true); // Number of direction changes
    Writer<size_t> w_CC(DataDir + "CC.dat", true); // Number of direction changes


    std::vector<size_t> nDC(nBB.size(), 0);

    // Main loop
    for (int k = 0; k < M; ++k) {
	if (10 * k % M == 0 && k != 0)
	    std::cout << "Simulated " << k << " initial conditions." << std::endl;

	Particle<T> P = Table.getInitialState(0, -1, TIME.back());
	//Particle<T> P0 = P;
	//P.Vel = Point<T>(cos(thetas[k]), sin(thetas[k]));
	//P.Pos = Point<T>(dx + 1e-7, d / (T)2);
	//P.Pos = Point<T>(1e-7, (d+dy+dy) / (T)2);
	T theta = atan(P.Vel.y / P.Vel.x);

	// Get trajectory - discrete
	Trajectory<T> Traj = Table.getTrajectory(P, N);
	std::vector<T> x = Traj.getLiftedPositions_x();
	std::vector<T> y = Traj.getLiftedPositions_y();
	std::vector<Point<T>> Vel = Traj.getVelocities();
	std::vector<T> t = Traj.getLiftedFreePaths();
	std::vector<Point<long int>> cid = Traj.getCellIDs();

	// curloop
	std::vector<int> FlowDirections = getFlowDirections(x, Vel, cid, dy, dx);

	for (size_t i = 0; i < nBB.size(); ++i) {
	    std::vector<Point<T>> lambda = getLaminarSteps(t, x, FlowDirections, minstep, nBB[i]);

	    nDC[i] += lambda.size();
	    if (lambda.size() > 0) {
		w_CC.WriteRowVector(ConsecutiveCounts(FlowDirections));
		//WriteColVector<T>(DataDir + "n" + std::to_string(nBB[i]) + "_" + std::to_string(k) + "--t.dat", t);
		//WriteColVector<T>(DataDir + "n" + std::to_string(nBB[i]) + "_" + std::to_string(k) + "--x.dat", x);
		//WriteColVector<T>(DataDir + "n" + std::to_string(nBB[i]) + "_" + std::to_string(k) + "--y.dat", y);
		
		std::vector<T> tt(lambda.size(), 0);
		std::vector<T> th(lambda.size(), 0);
		for (size_t j = 0; j < lambda.size(); ++j) {
		    th[j] = theta;
		    tt[j] = P.TotalTime;
		}
		std::vector<T> lx = getCoordFromVectorPoint(lambda, "x");
		std::vector<T> ly = getCoordFromVectorPoint(lambda, "y");

		if (nDC[i] < 100000)
		    w_JUMP[i].WriteColVectors(lx, ly);//, tt, th);
	    }
	}
    }
    w_DC.WriteRowVector(nDC);
}
