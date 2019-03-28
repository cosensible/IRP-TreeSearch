#include "Solver.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>

#include <cmath>

#include "CsvReader.h"
#include "MpSolver.h"
#include "CachedTspSolver.h"


using namespace std;


namespace szx {

#pragma region Solver::Cli
	int Solver::Cli::run(int argc, char * argv[]) {
		Log(LogSwitch::Szx::Cli) << "parse command line arguments." << endl;
		Set<String> switchSet;
		Map<String, char*> optionMap({ // use string as key to compare string contents instead of pointers.
			{ InstancePathOption(), nullptr },
			{ SolutionPathOption(), nullptr },
			{ RandSeedOption(), nullptr },
			{ TimeoutOption(), nullptr },
			{ MaxIterOption(), nullptr },
			{ JobNumOption(), nullptr },
			{ RunIdOption(), nullptr },
			{ EnvironmentPathOption(), nullptr },
			{ ConfigPathOption(), nullptr },
			{ LogPathOption(), nullptr }
			});

		for (int i = 1; i < argc; ++i) { // skip executable name.
			auto mapIter = optionMap.find(argv[i]);
			if (mapIter != optionMap.end()) { // option argument.
				mapIter->second = argv[++i];
			}
			else { // switch argument.
				switchSet.insert(argv[i]);
			}
		}

		Log(LogSwitch::Szx::Cli) << "execute commands." << endl;
		if (switchSet.find(HelpSwitch()) != switchSet.end()) {
			cout << HelpInfo() << endl;
		}

		if (switchSet.find(AuthorNameSwitch()) != switchSet.end()) {
			cout << AuthorName() << endl;
		}

		Solver::Environment env;
		env.load(optionMap);
		if (env.instPath.empty() || env.slnPath.empty()) { return -1; }

		Solver::Configuration cfg;
		cfg.load(env.cfgPath);

		Log(LogSwitch::Szx::Input) << "load instance " << env.instPath << " (seed=" << env.randSeed << ")." << endl;
		Problem::Input input;
		if (!input.load(env.instPath)) { return -1; }

		Solver solver(input, env, cfg);
		solver.solve();

		pb::Submission submission;
		submission.set_thread(to_string(env.jobNum));
		submission.set_instance(env.friendlyInstName());
		submission.set_duration(to_string(solver.timer.elapsedSeconds()) + "s");
		submission.set_obj(solver.output.totalCost);

		solver.output.save(env.slnPath, submission);
#if SZX_DEBUG
		solver.output.save(env.solutionPathWithTime(), submission);
		solver.record();
#endif // SZX_DEBUG

		return 0;
	}
#pragma endregion Solver::Cli

#pragma region Solver::Environment
	void Solver::Environment::load(const Map<String, char*> &optionMap) {
		char *str;

		str = optionMap.at(Cli::EnvironmentPathOption());
		if (str != nullptr) { loadWithoutCalibrate(str); }

		str = optionMap.at(Cli::InstancePathOption());
		if (str != nullptr) { instPath = str; }

		str = optionMap.at(Cli::SolutionPathOption());
		if (str != nullptr) { slnPath = str; }

		str = optionMap.at(Cli::RandSeedOption());
		if (str != nullptr) { randSeed = atoi(str); }

		str = optionMap.at(Cli::TimeoutOption());
		if (str != nullptr) { msTimeout = static_cast<Duration>(atof(str) * Timer::MillisecondsPerSecond); }

		str = optionMap.at(Cli::MaxIterOption());
		if (str != nullptr) { maxIter = atoi(str); }

		str = optionMap.at(Cli::JobNumOption());
		if (str != nullptr) { jobNum = atoi(str); }

		str = optionMap.at(Cli::RunIdOption());
		if (str != nullptr) { rid = str; }

		str = optionMap.at(Cli::ConfigPathOption());
		if (str != nullptr) { cfgPath = str; }

		str = optionMap.at(Cli::LogPathOption());
		if (str != nullptr) { logPath = str; }

		calibrate();
	}

	void Solver::Environment::load(const String &filePath) {
		loadWithoutCalibrate(filePath);
		calibrate();
	}

	void Solver::Environment::loadWithoutCalibrate(const String &filePath) {
		// EXTEND[szx][8]: load environment from file.
		// EXTEND[szx][8]: check file existence first.
	}

	void Solver::Environment::save(const String &filePath) const {
		// EXTEND[szx][8]: save environment to file.
	}
	void Solver::Environment::calibrate() {
		// adjust thread number.
		int threadNum = thread::hardware_concurrency();
		if ((jobNum <= 0) || (jobNum > threadNum)) { jobNum = threadNum; }

		// adjust timeout.
		msTimeout -= Environment::SaveSolutionTimeInMillisecond;
	}
#pragma endregion Solver::Environment

#pragma region Solver::Configuration
	void Solver::Configuration::load(const String &filePath) {
		// EXTEND[szx][5]: load configuration from file.
		// EXTEND[szx][8]: check file existence first.
	}

	void Solver::Configuration::save(const String &filePath) const {
		// EXTEND[szx][5]: save configuration to file.
	}
#pragma endregion Solver::Configuration

#pragma region Solver
	bool Solver::solve() {
		init();

		int workerNum = (max)(1, env.jobNum / cfg.threadNumPerWorker);
		cfg.threadNumPerWorker = env.jobNum / workerNum;
		List<Solution> solutions(workerNum, Solution(this));
		List<bool> success(workerNum);

		Log(LogSwitch::Szx::Framework) << "launch " << workerNum << " workers." << endl;
		List<thread> threadList;
		threadList.reserve(workerNum);
		for (int i = 0; i < workerNum; ++i) {
			// TODO[szx][2]: as *this is captured by ref, the solver should support concurrency itself, i.e., data members should be read-only or independent for each worker.
			// OPTIMIZE[szx][3]: add a list to specify a series of algorithm to be used by each threads in sequence.
			threadList.emplace_back([&, i]() { success[i] = optimize(solutions[i], i); });
		}
		for (int i = 0; i < workerNum; ++i) { threadList.at(i).join(); }

		Log(LogSwitch::Szx::Framework) << "collect best result among all workers." << endl;
		int bestIndex = -1;
		double bestValue = 0;
		for (int i = 0; i < workerNum; ++i) {
			if (!success[i]) { continue; }
			Log(LogSwitch::Szx::Framework) << "worker " << i << " got " << solutions[i].totalCost << endl;
			if (solutions[i].totalCost <= bestValue) { continue; }
			bestIndex = i;
			bestValue = solutions[i].totalCost;
		}

		env.rid = to_string(bestIndex);
		if (bestIndex < 0) { return false; }
		output = solutions[bestIndex];
		return true;
	}

	void Solver::record() const {
#if SZX_DEBUG
		int generation = 0;

		ostringstream log;

		System::MemoryUsage mu = System::peakMemoryUsage();

		// load reference results.
		CsvReader cr;
		ifstream ifs(Environment::DefaultInstanceDir() + "Baseline.csv");
		if (!ifs.is_open()) { return; }
		const List<CsvReader::Row> &rows(cr.scan(ifs));
		ifs.close();
		String bestObj, refObj, refTime;
		for (auto r = rows.begin(); r != rows.end(); ++r) {
			if (env.friendlyInstName() != r->front()) { continue; }
			bestObj = (*r)[1];
			refObj = (*r)[2];
			refTime = (*r)[3];
			//double opt = stod(bestObj);
			break;
		}

		double checkerObj = -1;
		bool feasible = check(checkerObj);
		double objDiff = round(output.totalCost * Problem::CheckerObjScale - checkerObj) / Problem::CheckerObjScale;

		// record basic information.
		log << env.friendlyLocalTime() << ","
			<< env.rid << ","
			<< env.instPath << ","
			<< feasible << "," << objDiff << ","
			<< output.totalCost << ","
			<< bestObj << ","
			<< refObj << ","
			<< timer.elapsedSeconds() << ","
			<< refTime << ","
			<< mu.physicalMemory << "," << mu.virtualMemory << ","
			<< env.randSeed << ","
			<< cfg.toBriefStr() << ","
			<< generation << "," << iteration;

		// record solution vector.
		// EXTEND[szx][2]: save solution in log.
		log << endl;

		// append all text atomically.
		static mutex logFileMutex;
		lock_guard<mutex> logFileGuard(logFileMutex);

		ofstream logFile(env.logPath, ios::app);
		logFile.seekp(0, ios::end);
		if (logFile.tellp() <= 0) {
			logFile << "Time,ID,Instance,Feasible,ObjMatch,Cost,MinCost,RefCost,Duration,RefDuration,PhysMem,VirtMem,RandSeed,Config,Generation,Iteration,Solution" << endl;
		}
		logFile << log.str();
		logFile.close();
#endif // SZX_DEBUG
	}

	bool Solver::check(double &checkerObj) const {
#if SZX_DEBUG
		enum CheckerFlag {
			IoError = 0x0,
			FormatError = 0x1,
			MultipleVisitsError = 0x2,
			UnmatchedLoadDeliveryError = 0x4,
			ExceedCapacityError = 0x8,
			RunOutOfStockError = 0x10
		};

		int errorCode = System::exec("Checker.exe " + env.instPath + " " + env.solutionPathWithTime());
		if (errorCode > 0) {
			checkerObj = errorCode;
			return true;
		}
		errorCode = ~errorCode;
		if (errorCode == CheckerFlag::IoError) { Log(LogSwitch::Checker) << "IoError." << endl; }
		if (errorCode & CheckerFlag::FormatError) { Log(LogSwitch::Checker) << "FormatError." << endl; }
		if (errorCode & CheckerFlag::MultipleVisitsError) { Log(LogSwitch::Checker) << "MultipleVisitsError." << endl; }
		if (errorCode & CheckerFlag::UnmatchedLoadDeliveryError) { Log(LogSwitch::Checker) << "UnmatchedLoadDeliveryError." << endl; }
		if (errorCode & CheckerFlag::ExceedCapacityError) { Log(LogSwitch::Checker) << "ExceedCapacityError." << endl; }
		if (errorCode & CheckerFlag::RunOutOfStockError) { Log(LogSwitch::Checker) << "RunOutOfStockError." << endl; }
		return false;
#else
		checkerObj = 0;
		return true;
#endif // SZX_DEBUG
	}

	void Solver::init() {
		nodeNum = input.nodes_size(), periodNum = input.periodnum();
		aux.routingCost.init(nodeNum, nodeNum);
		aux.routingCost.reset(Arr2D<Price>::ResetOption::AllBits0);
		aux.visits.init(periodNum, nodeNum);
		aux.visits.reset(Arr2D<ID>::ResetOption::AllBits0);
		aux.tours.init(periodNum);
		H1.resize(BitSize); H2.resize(BitSize); H3.resize(BitSize);

		ID n = 0;
		for (auto i = input.nodes().begin(); i != input.nodes().end(); ++i, ++n) {
			ID m = 0;
			for (auto j = input.nodes().begin(); j != i; ++j, ++m) {
				double value = round(hypot(i->x() - j->x(), i->y() - j->y()));
				aux.routingCost[n][m] = aux.routingCost[m][n] = value;
			}
		}

		aux.initHoldingCost = 0;
		for (auto i = input.nodes().begin(); i != input.nodes().end(); ++i) {
			aux.initHoldingCost += i->holdingcost() * i->initquantity();
		}
	}

	ID Solver::hash1(const Arr2D<ID> &visits) {
		unsigned long long sum = 0;
		for (ID p = 0; p < periodNum; ++p) {
			for (ID n = 0; n < nodeNum; ++n) {
				if (visits[p][n]) {
					sum += static_cast<int>(std::pow(p*nodeNum + n, gamma1));
				}
			}
		}
		return sum % BitSize;
	}

	ID Solver::hash2(const Arr2D<ID> &visits) {
		unsigned long long sum = 0;
		for (ID p = 0; p < periodNum; ++p) {
			for (ID n = 0; n < nodeNum; ++n) {
				if (visits[p][n]) {
					sum += static_cast<int>(std::pow(p*nodeNum + n, gamma2));
				}
			}
		}
		return sum % BitSize;
	}

	ID Solver::hash3(const Arr2D<ID> &visits) {
		unsigned long long sum = 0;
		for (ID p = 0; p < periodNum; ++p) {
			for (ID n = 0; n < nodeNum; ++n) {
				if (visits[p][n]) {
					sum += static_cast<int>(std::pow(p*nodeNum + n, gamma3));
				}
			}
		}
		return sum % BitSize;
	}

	bool Solver::isTabu(unsigned hv1, unsigned hv2, unsigned hv3, ID to1, ID to0) {
		hv1 += (static_cast<unsigned>(std::pow(to1, gamma1)) - static_cast<unsigned>(std::pow(to0, gamma1)));
		hv2 += (static_cast<unsigned>(std::pow(to1, gamma2)) - static_cast<unsigned>(std::pow(to0, gamma2)));
		hv3 += (static_cast<unsigned>(std::pow(to1, gamma3)) - static_cast<unsigned>(std::pow(to0, gamma3)));
		return (H1[hv1] && H2[hv2] && H3[hv3]);
	}

	void Solver::execTabu(ID to1, ID to0) {
		hashValue1 += (static_cast<unsigned>(std::pow(to1, gamma1)) - static_cast<unsigned>(std::pow(to0, gamma1)));
		hashValue2 += (static_cast<unsigned>(std::pow(to1, gamma2)) - static_cast<unsigned>(std::pow(to0, gamma2)));
		hashValue3 += (static_cast<unsigned>(std::pow(to1, gamma3)) - static_cast<unsigned>(std::pow(to0, gamma3)));
		H1[hashValue1] = H2[hashValue2] = H3[hashValue3] = 1;
	}

	void Solver::sln2Visits(Solution &sln) {
		aux.bestCost = sln.totalCost;
		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < input.vehicles_size(); ++v) {
				auto &delivs(*sln.mutable_periodroutes(p)->mutable_vehicleroutes(v)->mutable_deliveries());
				for (auto n = delivs.cbegin(); n != delivs.cend(); ++n) {
					aux.visits[p][n->node()] = 1;
				}
			}
		}
	}

	int Solver::buildDelNeigh(Solution &sln, Arr2D<ID> &visits) {
		aux.delNeigh.clear();
		for (ID p = 0; p < periodNum; ++p) {
			for (ID n = 1; n < nodeNum; ++n) {
				if (visits[p][n]) {
					visits[p][n] = 0;
					Price totalCost = callModel(sln, visits);
					if (totalCost >= 0 && totalCost <= aux.bestCost) {
						aux.delNeigh.push_back({ totalCost,p*nodeNum + n });
					}
					visits[p][n] = 1;
				}
			}
		}
		sort(aux.delNeigh.begin(), aux.delNeigh.end());	// 选最优的邻居解
		return aux.delNeigh.size();
	}

	int Solver::buildAddNeigh(Solution &sln, Arr2D<ID> &visits) {
		aux.addNeigh.clear();
		const auto &nodes(*input.mutable_nodes());
		for (ID p = 0; p < periodNum; ++p) {
			for (ID n = 1; n < nodeNum; ++n) {
				if (!visits[p][n] && nodes[n].holdingcost() < nodes[0].holdingcost()) {
					visits[p][n] = visits[p][0] = 1;
					Price totalCost = callModel(sln, visits);
					if (totalCost >= 0) {
						aux.addNeigh.push_back({ totalCost,p*nodeNum + n });
					}
					visits[p][n] = 0;
				}
			}
		}
		sort(aux.addNeigh.begin(), aux.addNeigh.end());	// 选最优的邻居解
		return aux.addNeigh.size();
	}

	int Solver::buildSwapNeigh(Solution &sln, Arr2D<ID> &visits) {
		aux.swapNeigh.clear();
		//for (ID p = 0; p < periodNum; ++p) {
		//	List<ID> N0, N1;
		//	for (ID n = 1; n < nodeNum; ++n) {
		//		if (visits[p][n]) { N1.push_back(n); }
		//		else { N0.push_back(n); }
		//	}
		//	for (ID &n0 : N0) {
		//		visits[p][n0] = 1;
		//		for (ID &n1 : N1) {
		//			visits[p][n1] = 0;
		//			Price totalCost = callModel(sln, visits);
		//			if (totalCost >= 0 && totalCost <= aux.bestCost) {
		//				aux.swapNeigh.push_back(SwapActor(p*nodeNum + n0, p*nodeNum + n1, totalCost));
		//			}
		//			visits[p][n1] = 1;
		//		}
		//		visits[p][n0] = 0;
		//	}
		//}

		for (ID n = 1; n < nodeNum; ++n) {
			List<ID> P0, P1;
			for (ID p = 0; p < periodNum; ++p) {
				if (visits[p][n]) { P1.push_back(p); }
				else { P0.push_back(p); }
			}
			for (ID &p0 : P0) {
				visits[p0][n] = visits[p0][0] = 1;	// 若该周期仓库为0，则仓库也要变为1
				for (ID &p1 : P1) {
					visits[p1][n] = 0;
					if (!isTabu(hashValue1,hashValue2,hashValue3, p0*nodeNum + n, p1*nodeNum + n)) {
						Price totalCost = callModel(sln, visits);
						if (totalCost >= 0 && totalCost <= aux.bestCost) {
							aux.swapNeigh.push_back({ totalCost,p0*nodeNum + n,p1*nodeNum + n });
						}
					}
					visits[p1][n] = 1;
				}
				visits[p0][n] = 0;
			}
		}
		sort(aux.swapNeigh.begin(), aux.swapNeigh.end());	// 选最优的邻居解
		return aux.swapNeigh.size();
	}

	int Solver::buildTabuNeigh(Solution &sln, Arr2D<ID> &visits) {
		aux.tabuNeigh.clear();
		//for (ID p = 0; p < periodNum; ++p) {
		//	List<ID> N0, N1;
		//	for (ID n = 1; n < nodeNum; ++n) {
		//		if (visits[p][n]) { N1.push_back(n); }
		//		else { N0.push_back(n); }
		//	}
		//	for (ID &n0 : N0) {
		//		visits[p][n0] = 1;
		//		for (ID &n1 : N1) {
		//			visits[p][n1] = 0;
		//			if (!isTabu(hashValue1, hashValue2, hashValue3, p*nodeNum + n0, p*nodeNum + n1)) {
		//				Price totalCost = callModel(sln, visits);
		//				if (totalCost >= 0) {	// 模型有解
		//					if (0 == aux.tabuNeigh.size() || totalCost == std::get<0>(aux.tabuNeigh[0])) {
		//						aux.tabuNeigh.push_back({ totalCost,p*nodeNum + n0,p*nodeNum + n1 });
		//					}
		//					else if (totalCost < std::get<0>(aux.tabuNeigh[0])) {
		//						aux.tabuNeigh.clear();
		//						aux.tabuNeigh.push_back({ totalCost,p*nodeNum + n0,p*nodeNum + n1 });
		//					}
		//				}
		//			}
		//			visits[p][n1] = 1;
		//		}
		//		visits[p][n0] = 0;
		//	}
		//}

		for (ID n = 1; n < nodeNum; ++n) {
			List<ID> P0, P1;
			for (ID p = 0; p < periodNum; ++p) {
				if (visits[p][n]) { P1.push_back(p); }
				else { P0.push_back(p); }
			}
			for (ID &p0 : P0) {
				visits[p0][n] = visits[p0][0] = 1;	// 若该周期仓库为0，则仓库也要变为1
				for (ID &p1 : P1) {
					visits[p1][n] = 0;
					if (!isTabu(hashValue1, hashValue2, hashValue3, p0*nodeNum + n, p1*nodeNum + n)) {
						Price totalCost = callModel(sln, visits);
						if (totalCost >= 0) {	// 模型有解
							if (0 == aux.tabuNeigh.size() || totalCost == std::get<0>(aux.tabuNeigh[0])) {
								aux.tabuNeigh.push_back({ totalCost,p0*nodeNum + n,p1*nodeNum + n });
							}
							else if (totalCost < std::get<0>(aux.tabuNeigh[0])) {
								aux.tabuNeigh.clear();
								aux.tabuNeigh.push_back({ totalCost,p0*nodeNum + n,p1*nodeNum + n });
							}
						}
					}
					visits[p1][n] = 1;
				}
				visits[p0][n] = 0;
			}
		}
		return aux.tabuNeigh.size();
	}

	void Solver::swapTabuSearch(Solution &sln, Arr2D<ID> &visits) {
		Log(LogSwitch::Szx::STS) << "Tabu Search" << endl;
		int tabuNeighSize = 0, step = 0;
		Price cost; ID to1, to0;
		Arr2D<ID> tabuVisits(visits);
		hashValue1 = hash1(tabuVisits), hashValue2 = hash2(tabuVisits), hashValue3 = hash3(tabuVisits);
		H1[hashValue1] = H2[hashValue2] = H3[hashValue3] = 1;	// 禁忌初始解
		while ((step++) < alpha && (tabuNeighSize = buildTabuNeigh(sln, tabuVisits))) {
			std::tie(cost, to1, to0) = aux.tabuNeigh[rand.pick(tabuNeighSize)];
			tabuVisits[to1 / nodeNum][to1 % nodeNum] = 1;
			tabuVisits[to0 / nodeNum][to0 % nodeNum] = 0;
			execTabu(to1, to0);	// 禁忌该解
			if (cost < aux.bestCost) {
				step = 0;
				aux.bestCost = cost;
				aux.visits = tabuVisits;
				Log(LogSwitch::Szx::STS) << "By STS, opt=" << aux.bestCost << endl;
			}
			if (localSearch(sln, tabuVisits)) {
				step = 0;
				aux.visits = tabuVisits;
			}
		}
	}

	void Solver::finalSearch(Solution &sln) {
		Log(LogSwitch::Szx::STS) << "Final Search" << endl;
		Price cost; ID to1;
		Arr2D<ID> finalVisits(aux.visits);
		buildAddNeigh(sln, finalVisits);
		ID maxAddNum = aux.addNeigh.size() > 3 ? 3 : aux.addNeigh.size();
		for (ID step = 0; step < maxAddNum; ++step) {
			std::tie(cost, to1) = aux.addNeigh[step];
			finalVisits[to1 / nodeNum][to1 % nodeNum] = 1;
			if (cost < aux.bestCost) {
				aux.bestCost = cost;
				aux.visits = finalVisits;
				Log(LogSwitch::Szx::STS) << "By FS, opt=" << aux.bestCost << endl;
			}
			swapTabuSearch(sln, finalVisits);
			finalVisits[to1 / nodeNum][to1 % nodeNum] = 0;
		}
	}

	bool Solver::localSearch(Solution &sln, Arr2D<ID> &visits) {
		int delNeighSize = 0, swapNeighSize = 0;
		bool isImproved = false;
		do {
			while (delNeighSize = buildDelNeigh(sln, visits)) {
				Log(LogSwitch::Szx::LS) << "delNeighSize = " << delNeighSize << endl;
				//ID to0; std::tie(aux.bestCost, to0) = aux.delNeigh[rand.pick(delNeighSize)]; // 在邻域随机选择
				ID to0; std::tie(aux.bestCost, to0) = aux.delNeigh[0];	// 选择最佳邻居解
				visits[to0 / nodeNum][to0 % nodeNum] = 0;
				isImproved = true;
				Log(LogSwitch::Szx::LS) << "By del, opt=" << aux.bestCost << endl;
			}
			while (swapNeighSize = buildSwapNeigh(sln, visits)) {
				Log(LogSwitch::Szx::LS) << "swapNeighSize = " << swapNeighSize << endl;
				//ID to1, to0; std::tie(aux.bestCost, to1, to0) = aux.swapNeigh[rand.pick(swapNeighSize)]; // 在邻域随机选择
				ID to1, to0; std::tie(aux.bestCost, to1, to0) = aux.swapNeigh[0];	// 选择最佳邻居解
				visits[to1 / nodeNum][to1 % nodeNum] = 1;
				visits[to0 / nodeNum][to0 % nodeNum] = 0;
				execTabu(to1, to0);	// 禁忌该解
				isImproved = true;
				Log(LogSwitch::Szx::LS) << "By swap, opt=" << aux.bestCost << endl;
			}
		} while (delNeighSize || swapNeighSize);
		return isImproved;
	}

	void Solver::getBestSln(Solution &sln) {
		callModel(sln, aux.visits, true);
	}

	bool Solver::optimize(Solution &sln, ID workerId) {
		Log(LogSwitch::Szx::Framework) << "worker " << workerId << " starts." << endl;

		sln.init(periodNum, input.vehicles_size(), Problem::MaxCost);
		iteratedModel(sln);

		localSearch(sln, aux.visits);

		swapTabuSearch(sln,aux.visits);

		finalSearch(sln);

		getBestSln(sln);

		Log(LogSwitch::Szx::Framework) << "worker " << workerId << " ends." << endl;
		return true;
	}

	void Solver::iteratedModel(Solution &sln) {
		ID vehicleNum = input.vehicles_size();
		const auto &vehicles(*input.mutable_vehicles());
		const auto &nodes(*input.mutable_nodes());

		MpSolver::Configuration mpCfg(MpSolver::InternalSolver::GurobiMip, env.timeoutInSecond(), true, false);
		MpSolver mp(mpCfg); mp.setMaxThread(4);

		// delivery[p, v, n] is the quantity delivered to node n at period p by vehicle v.
		Arr2D<Arr<Dvar>> delivery(periodNum, vehicleNum, Arr<Dvar>(nodeNum));
		// x[p, v, n, m] is true if the edge from node n to node m is visited at period p by vehicle v.
		Arr2D<Arr2D<Dvar>> x(periodNum, vehicleNum, Arr2D<Dvar>(nodeNum, nodeNum));

		// quantityLevel[n, p] is the rest quantity of node n at period p after the delivery and consumption have happened.
		Arr2D<Expr> quantityLevel(nodeNum, periodNum);
		// degrees[n] is the sum of in-coming edges.
		Arr2D<Arr<Expr>> degrees(periodNum, vehicleNum, Arr<Expr>(nodeNum));

		// add decision variables.
		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				for (ID n = 0; n < input.depotnum(); ++n) {
					Quantity capacity = min(input.vehicles(v).capacity(), input.nodes(n).capacity());
					delivery[p][v][n] = mp.addVar(MpSolver::VariableType::Real, -capacity, 0);
				}
				for (ID n = input.depotnum(); n < nodeNum; ++n) {
					Quantity capacity = min(input.vehicles(v).capacity(), input.nodes(n).capacity());
					delivery[p][v][n] = mp.addVar(MpSolver::VariableType::Real, 0, capacity);
				}
				Arr2D<Dvar> &xpv(x.at(p, v));
				for (ID n = 0; n < nodeNum; ++n) {
					for (ID m = 0; m < nodeNum; ++m) {
						if (n == m) { continue; }
						xpv.at(n, m) = mp.addVar(MpSolver::VariableType::Bool, 0, 1);
					}
				}
			}
		}

		// add constraints.
		for (ID n = 0; n < nodeNum; ++n) {
			const auto &node(input.nodes(n));
			Expr quantity = node.initquantity();
			for (ID p = 0; p < periodNum; ++p) {
				for (ID v = 0; v < vehicleNum; ++v) {
					quantity += delivery[p][v][n];
				}
				// node capacity constraint.
				mp.addConstraint(quantity <= node.capacity());
				quantity -= node.demands(p);
				mp.addConstraint(0 <= quantity);
				quantityLevel[n][p] = quantity;
			}
		}

		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				Expr quantity;
				for (ID n = 0; n < nodeNum; ++n) {
					quantity += delivery[p][v][n];
				}
				// quantity matching constraint.
				mp.addConstraint(quantity == 0);
			}
		}

		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				Arr2D<Dvar> &xpv(x.at(p, v));
				for (ID n = 0; n < nodeNum; ++n) {
					Expr inDegree;
					for (ID m = 0; m < n; ++m) { inDegree += xpv.at(m, n); }
					for (ID m = n + 1; m < nodeNum; ++m) { inDegree += xpv.at(m, n); }
					degrees[p][v][n] = inDegree;
					Expr outDegree;
					for (ID m = 0; m < n; ++m) { outDegree += xpv.at(n, m); }
					for (ID m = n + 1; m < nodeNum; ++m) { outDegree += xpv.at(n, m); }
					// path connectivity constraint.
					mp.addConstraint(inDegree == outDegree); // OPTIMIZE[szx][0]: use undirected graph version? (degree == 2)
					// delivery precondition constraint.
					Quantity capacity = min(input.vehicles(v).capacity(), input.nodes(n).capacity());
					double quantityCoef = (n >= input.depotnum()) ? 1 : -1;
					mp.addConstraint(quantityCoef * delivery[p][v][n] <= capacity * inDegree);
					if (n >= input.depotnum()) {
						// visit precondition constraint.
						mp.addConstraint(delivery[p][v][n] >= inDegree); // OPTIMIZE[szx][2]: omit it since it will be satisfied automatically?
					}
					// maximal visit constraint.
					mp.addConstraint(inDegree <= 1); // OPTIMIZE[szx][2]: omit it since it will be satisfied automatically?
				}
			}
		}

		// add objective.
		Expr obj;
		Expr holdingCost = aux.initHoldingCost;
		for (ID n = 0; n < nodeNum; ++n) {
			const auto &node(input.nodes(n));
			for (ID p = 0; p < periodNum; ++p) {
				holdingCost += (node.holdingcost() * quantityLevel.at(n, p));
			}
		}
		Expr routingCost;
		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				Arr2D<Dvar> &xpv(x.at(p, v));
				for (ID n = 0; n < nodeNum; ++n) {
					for (ID m = 0; m < nodeNum; ++m) {
						if (n == m) { continue; }
						routingCost += (aux.routingCost.at(n, m) * xpv.at(n, m));
					}
				}
			}
		}
		obj = holdingCost + routingCost;
		mp.addObjective(obj, MpSolver::OptimaOrientation::Minimize, 0, 0, 0, env.timeoutInSecond());

		// add callbacks.
		auto subTourHandler = [&](MpSolver::MpEvent &e) {
			enum EliminationPolicy { // OPTIMIZE[szx][0]: first sub-tour, best sub-tour or all sub-tours?
				NoSubTour = 0x0,
				AllSubTours = 0x1,
				FirstSubTour = 0x2,
				BestSubTour = 0x4
			};
			EliminationPolicy policy = EliminationPolicy::BestSubTour;

			List<ID> bestTour; // tour with least nodes/hops.
			List<ID> tour;
			tour.reserve(nodeNum);
			Arr<bool> visited(nodeNum);
			for (ID p = 0; p < periodNum; ++p) {
				for (ID v = 0; v < vehicleNum; ++v) {
					Arr2D<Dvar> &xpv(x.at(p, v));
					tour.clear();
					visited.reset(Arr<bool>::ResetOption::AllBits0);
					for (ID s = 0; s < nodeNum; ++s) { // check if there is a route start from each node.
						if (visited[s]) { continue; }
						ID prev = s;
						do {
							for (ID n = 0; n < nodeNum; ++n) {
								if (prev == n) { continue; }
								if (!e.isTrue(xpv.at(prev, n))) { continue; }
								if (s >= input.depotnum()) { tour.push_back(n); } // the sub-tour containing depots should not be eliminated.
								prev = n;
								visited[n] = true;
								break;
							}
						} while (prev != s);
						if (tour.empty()) { continue; }

						if (policy & (EliminationPolicy::AllSubTours | EliminationPolicy::FirstSubTour)) {
							Expr edges;
							for (auto n = tour.begin(); n != tour.end(); prev = *n, ++n) {
								edges += xpv.at(prev, *n);
							}
							e.addLazy(edges <= static_cast<double>(tour.size() - 1));
							if (policy & EliminationPolicy::FirstSubTour) { break; }
						}

						if (bestTour.empty() || (tour.size() < bestTour.size())) { swap(bestTour, tour); }
					}
					if ((policy & EliminationPolicy::BestSubTour) && !bestTour.empty()) {
						Expr edges;
						ID prev = bestTour.back();
						for (auto n = bestTour.begin(); n != bestTour.end(); prev = *n, ++n) {
							edges += xpv.at(prev, *n);
						}
						e.addLazy(edges <= static_cast<double>(bestTour.size() - 1));
					}
				}
			}
		};

		static const String TspCacheDir("TspCache/");
		System::makeSureDirExist(TspCacheDir);
		CachedTspSolver tspSolver(nodeNum, TspCacheDir + env.friendlyInstName() + ".csv");

		Solution curSln; // current solution.
		curSln.init(periodNum, vehicleNum);
		auto nodeSetHandler = [&](MpSolver::MpEvent &e) {
			//if (e.getObj() > sln.totalCost) { e.stop(); } // there could be bad heuristic solutions.

			// OPTIMIZE[szx][0]: check the bound and only apply this to the optimal sln.

			lkh::CoordList2D coords; // OPTIMIZE[szx][3]: use adjacency matrix to avoid re-calculation and different rounding?
			coords.reserve(nodeNum);
			List<ID> nodeIdMap(nodeNum);
			List<bool> containNode(nodeNum);
			lkh::Tour tour;
			Expr nodeDiff;

			curSln.totalCost = 0;
			for (ID p = 0; p < periodNum; ++p) {
				auto &periodRoute(*curSln.mutable_periodroutes(p));
				for (ID v = 0; v < vehicleNum; ++v) {
					Arr2D<Dvar> &xpv(x.at(p, v));
					coords.clear();
					fill(containNode.begin(), containNode.end(), false);
					for (ID n = 0; n < nodeNum; ++n) {
						bool visited = false;
						for (ID m = 0; m < nodeNum; ++m) {
							if (n == m) { continue; }
							if (!e.isTrue(xpv.at(n, m))) { continue; }
							nodeIdMap[coords.size()] = n;
							containNode[n] = true;
							coords.push_back(lkh::Coord2D(nodes[n].x() * Precision, nodes[n].y() * Precision));
							visited = true;
							break;
						}
						nodeDiff += (visited ? (1 - degrees[p][v][n]) : degrees[p][v][n]);
					}
					auto &route(*periodRoute.mutable_vehicleroutes(v));
					route.clear_deliveries();
					if (coords.size() > 2) { // repair the relaxed solution.
						tspSolver.solve(tour, containNode, coords, [&](ID n) { return nodeIdMap[n]; });
					}
					else if (coords.size() == 2) { // trivial cases.
						tour.nodes.resize(2);
						tour.nodes[0] = nodeIdMap[0];
						tour.nodes[1] = nodeIdMap[1];
					}
					else {
						continue;
					}
					tour.nodes.push_back(tour.nodes.front());
					for (auto n = tour.nodes.begin(), m = n + 1; m != tour.nodes.end(); ++n, ++m) {
						auto &d(*route.add_deliveries());
						d.set_node(*m);
						d.set_quantity(lround(e.getValue(delivery[p][v][*m])));
						curSln.totalCost += aux.routingCost.at(*n, *m);
					}
				}
			}
			e.addLazy(nodeDiff >= 1);

			curSln.totalCost += e.getValue(holdingCost);
			if (curSln.totalCost < sln.totalCost) {
				Log(LogSwitch::Szx::Model) << "By nodeSetHandler, opt=" << curSln.totalCost << endl;
				std::swap(curSln, sln);
			}
		};

		//mp.setMipSlnEvent([&](MpSolver::MpEvent &e) {
		//	nodeSetHandler(e);
		//	subTourHandler(e);
		//});

		int iter = 0, mod = 3;
		mp.setMipSlnEvent([&](MpSolver::MpEvent &e) {
			if ((++iter) % mod) { subTourHandler(e); }
			else { nodeSetHandler(e); }
			
			//if ((++iter) % mod) { nodeSetHandler(e); }
			//else { subTourHandler(e); }
		});

		// solve.
		if (mp.optimize()) {
			curSln.totalCost = obj.getValue();
			for (ID p = 0; p < periodNum; ++p) {
				auto &periodRoute(*curSln.mutable_periodroutes(p));
				for (ID v = 0; v < vehicleNum; ++v) {
					Arr2D<Dvar> &xpv(x.at(p, v));
					auto &route(*periodRoute.mutable_vehicleroutes(v));
					route.clear_deliveries();
					for (ID s = 0; s < input.depotnum(); ++s) {
						ID prev = s;
						do {
							for (ID n = 0; n < nodeNum; ++n) {
								if (prev == n) { continue; }
								if (!mp.isTrue(xpv.at(prev, n))) { continue; }
								auto &d(*route.add_deliveries());
								d.set_node(n);
								d.set_quantity(lround(mp.getValue(delivery[p][v][n])));
								prev = n;
								break;
							}
						} while (prev != s);
					}
				}
			}
			if (curSln.totalCost < sln.totalCost) {
				Log(LogSwitch::Szx::Model) << "By subtourHandler, opt=" << curSln.totalCost << endl;
				std::swap(curSln, sln);
			}
		}

		sln2Visits(sln);	// 初始化 visits 和 bestCost
		// 禁忌初始解
		hashValue1 = hash1(aux.visits), hashValue2 = hash2(aux.visits), hashValue3 = hash3(aux.visits);
		H1[hashValue1] = H2[hashValue2] = H3[hashValue3] = 1;
	}

	// isBest 表示在记录最优解
	Price Solver::callModel(Solution &sln, Arr2D<int> &visits, bool isBest) {
		ID vehicleNum = input.vehicles_size();

		MpSolver::Configuration mpCfg(MpSolver::InternalSolver::GurobiMip, env.timeoutInSecond(), true, false);
		MpSolver mp(mpCfg); mp.setMaxThread(4);

		// delivery[p, v, n] is the quantity delivered to node n at period p by vehicle v.
		Arr2D<Arr<Dvar>> delivery(periodNum, vehicleNum, Arr<Dvar>(nodeNum));
		// quantityLevel[n, p] is the rest quantity of node n at period p after the delivery and consumption have happened.
		Arr2D<Expr> quantityLevel(nodeNum, periodNum);

		// add decision variables.
		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				for (ID n = 0; n < input.depotnum(); ++n) {
					Quantity capacity = min(input.vehicles(v).capacity(), input.nodes(n).capacity());
					delivery[p][v][n] = mp.addVar(MpSolver::VariableType::Real, -capacity, 0);
				}
				for (ID n = input.depotnum(); n < nodeNum; ++n) {
					Quantity capacity = min(input.vehicles(v).capacity(), input.nodes(n).capacity());
					delivery[p][v][n] = mp.addVar(MpSolver::VariableType::Real, 0, capacity);
				}
			}
		}

		// add constraints.
		for (ID n = 0; n < nodeNum; ++n) {
			const auto &node(input.nodes(n));
			Expr quantity = node.initquantity();
			for (ID p = 0; p < periodNum; ++p) {
				for (ID v = 0; v < vehicleNum; ++v) {
					quantity += delivery[p][v][n];
				}
				// node capacity constraint.
				mp.addConstraint(quantity <= node.capacity());
				quantity -= node.demands(p);
				mp.addConstraint(0 <= quantity);
				quantityLevel[n][p] = quantity;
			}
		}

		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				Expr quantity;
				for (ID n = 0; n < nodeNum; ++n) {
					quantity += delivery[p][v][n];

					Quantity capacity = min(input.vehicles(v).capacity(), input.nodes(n).capacity());
					double quantityCoef = (n >= input.depotnum()) ? 1 : -1;
					ostringstream name; name << "V" << p << n;
					mp.addConstraint(quantityCoef * delivery[p][v][n] <= capacity * visits[p][n], name.str());
				}
				// quantity matching constraint.
				mp.addConstraint(quantity == 0);
			}
		}

		// add objective.
		Expr holdingCost = aux.initHoldingCost;
		for (ID n = 0; n < nodeNum; ++n) {
			const auto &node(input.nodes(n));
			for (ID p = 0; p < periodNum; ++p) {
				holdingCost += (node.holdingcost() * quantityLevel.at(n, p));
			}
		}
		mp.addObjective(holdingCost, MpSolver::OptimaOrientation::Minimize, 0, 0, 0, env.timeoutInSecond());

		if (mp.optimize()) {
			Price totalCost = callLKH(visits, isBest) + mp.getObjectiveValue();
			if (!isBest) { return  totalCost; }
			// 不是记录最优解，直接返回，否则还要记录最优解
			sln.totalCost = totalCost;
			for (ID p = 0; p < periodNum; ++p) {
				for (ID v = 0; v < vehicleNum; ++v) {
					auto &route(*sln.mutable_periodroutes(p)->mutable_vehicleroutes(v));
					route.clear_deliveries();
					for (auto n = aux.tours[p].begin(); n != aux.tours[p].end(); ++n) {
						auto &d(*route.add_deliveries());
						d.set_node(*n);
						d.set_quantity(lround(mp.getValue(delivery[p][v][*n])));
					}
				}
			}
			return totalCost;
		}
		return -1;
	}

	Price Solver::callLKH(const Arr2D<ID> &visits, bool isBest) {
		const auto &nodes(*input.mutable_nodes());

		static const String TspCacheDir("TspCache/");
		System::makeSureDirExist(TspCacheDir);
		CachedTspSolver tspSolver(nodeNum, TspCacheDir + env.friendlyInstName() + ".csv");

		lkh::CoordList2D coords;
		coords.reserve(nodeNum);
		List<ID> nodeIdMap(nodeNum);
		List<bool> containNode(nodeNum);
		lkh::Tour tour;

		Price tourCost = 0.0;
		for (ID p = 0; p < periodNum; ++p) {
			coords.clear();
			fill(containNode.begin(), containNode.end(), false);
			for (ID n = 0; n < nodeNum; ++n) {
				if (visits[p][n]) {
					nodeIdMap[coords.size()] = n;
					containNode[n] = true;
					coords.push_back(lkh::Coord2D(nodes[n].x() * Precision, nodes[n].y() * Precision));
				}
			}
			if (coords.size() > 2) { // repair the relaxed solution.
				tspSolver.solve(tour, containNode, coords, [&](ID n) { return nodeIdMap[n]; });
			}
			else if (coords.size() == 2) { // trivial cases.
				tour.nodes.resize(2);
				tour.nodes[0] = nodeIdMap[0];
				tour.nodes[1] = nodeIdMap[1];
			}
			else {
				continue;
			}
			tour.nodes.push_back(tour.nodes.front());
			for (auto n = tour.nodes.begin(), m = n + 1; m != tour.nodes.end(); ++n, ++m) {
				tourCost += aux.routingCost.at(*n, *m);
				if (isBest) { aux.tours[p].push_back(*m); }
			}
		}
		return tourCost;
	}

#pragma endregion Solver

}
