#include "Solver.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>
#include <stack>
#include <cmath>
#include <queue>

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
		aux.routingCost.init(input.nodes_size(), input.nodes_size());
		aux.routingCost.reset(Arr2D<Price>::ResetOption::AllBits0);
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

	Price Solver::callModel(const Arr2D<int> &visits) {
		ID nodeNum = input.nodes_size();
		ID vehicleNum = input.vehicles_size();
		ID periodNum = input.periodnum();

		MpSolver::Configuration mpCfg(MpSolver::InternalSolver::GurobiMip, env.timeoutInSecond(), true, false);
		MpSolver mp(mpCfg);

		// delivery[p, v, n] is the quantity delivered to node n at period p by vehicle v.
		Arr2D<Arr<Dvar>> delivery(input.periodnum(), vehicleNum, Arr<Dvar>(nodeNum));
		// quantityLevel[n, p] is the rest quantity of node n at period p after the delivery and consumption have happened.
		Arr2D<Expr> quantityLevel(nodeNum, input.periodnum());

		// add decision variables.
		for (ID p = 0; p < input.periodnum(); ++p) {
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
			for (ID p = 0; p < input.periodnum(); ++p) {
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

		for (ID p = 0; p < input.periodnum(); ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				Expr quantity;
				for (ID n = 0; n < nodeNum; ++n) {
					quantity += delivery[p][v][n];

					Quantity capacity = min(input.vehicles(v).capacity(), input.nodes(n).capacity());
					double quantityCoef = (n >= input.depotnum()) ? 1 : -1;
					mp.addConstraint(quantityCoef * delivery[p][v][n] <= capacity * visits[p][n]);
				}
				// quantity matching constraint.
				mp.addConstraint(quantity == 0);
			}
		}
		// add objective.
		Expr holdingCost = aux.initHoldingCost;
		for (ID n = 0; n < nodeNum; ++n) {
			const auto &node(input.nodes(n));
			for (ID p = 0; p < input.periodnum(); ++p) {
				holdingCost += (node.holdingcost() * quantityLevel.at(n, p));
			}
		}
		mp.addObjective(holdingCost, MpSolver::OptimaOrientation::Minimize, 0, 0, 0, env.timeoutInSecond());
		if (mp.optimize()) {
			return mp.getObjectiveValue();
		}
		return -1;
	}

	Price Solver::callLKH(const Solution &sln, const Arr2D<ID> &visits, int p) {
		ID nodeNum = input.nodes_size();
		const auto &nodes(*input.mutable_nodes());

		static const String TspCacheDir("TspCache/");
		System::makeSureDirExist(TspCacheDir);
		CachedTspSolver tspSolver(nodeNum, TspCacheDir + env.friendlyInstName() + ".csv");

		lkh::CoordList2D coords;
		coords.reserve(nodeNum);
		List<ID> nodeIdMap(nodeNum);
		List<bool> containNode(nodeNum, false);
		lkh::Tour tour;

		for (ID n = 0; n < nodeNum; ++n) {
			if (visits[p][n]) {
				nodeIdMap[coords.size()] = n;
				containNode[n] = true;
				coords.push_back(lkh::Coord2D(nodes[n].x() * Precision, nodes[n].y() * Precision));
			}
		}
		if (coords.size() < 2) { return 0.0; }
		if (coords.size() > 2) { // repair the relaxed solution.
			tspSolver.solve(tour, containNode, coords, [&](ID n) { return nodeIdMap[n]; });
		}
		else { // trivial cases.
			tour.nodes.resize(2);
			tour.nodes[0] = nodeIdMap[0];
			tour.nodes[1] = nodeIdMap[1];
		}
		Price tourCost = 0.0;
		tour.nodes.push_back(tour.nodes.front());
		for (auto n = tour.nodes.begin(), m = n + 1; m != tour.nodes.end(); ++n, ++m) {
			tourCost += aux.routingCost.at(*n, *m);
		}
		return tourCost;
	}

	void Solver::treeSearch(Solution &sln, int depth, int iterTime) {
		ID nodeNum = input.nodes_size();
		ID vehicleNum = input.vehicles_size();
		ID periodNum = input.periodnum();
		const auto &nodes(*input.mutable_nodes());

		Arr2D<ID> visits(periodNum, nodeNum, 0);
		Arr<Price> toursCost(periodNum, 0.0);
		for (ID p = 0; p < periodNum; ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				const auto &delivs(*sln.mutable_periodroutes(p)->mutable_vehicleroutes(v)->mutable_deliveries());
				if (delivs.size() < 2) { continue; }
				ID s = delivs.rbegin()->node(), c = delivs.cbegin()->node();
				cout << s << "->" << c << endl;
				toursCost[p] += aux.routingCost.at(s, c);
				visits[p][c] = 1;
				for (auto n = delivs.cbegin(), m = n + 1; m != delivs.cend(); ++n, ++m) {
					toursCost[p] += aux.routingCost.at(n->node(), m->node());
					visits[p][m->node()] = 1;
				}
			}
		}

		//for (ID p = 0; p < periodNum; ++p) {
		//	for (ID n = 0; n < nodeNum; ++n) {
		//		cout << visits[p][n] << " ";
		//	}
		//	cout << endl;
		//}

		// Price 成本 vector<ID>：为树搜索路径，倒数第二个元素为下一轮迭代起始节点，最后一个元素为该叶节点的迭代周期
		queue<pair<Price, vector<ID>>> leaves;
		leaves.push({ 0,{-1,-1,0} });
		Price bestCost = sln.totalCost;
		vector<ID> bestChange;	// 记录最优解的变化路径
		int visitNum = periodNum * nodeNum, leafNum = 3;
		for (int step = 0; step < iterTime; ++step) {
			vector<pair<Price, vector<ID>>> topkLeaves;
			while (!leaves.empty()) {
				pair<Price, vector<ID>> leaf = leaves.front();
				if (leaf.second.back() > step) { break; }
				leaves.pop(); leaf.second.pop_back();
				int root = leaf.second.back();	// 下一轮迭代起点
				leaf.second.pop_back(); leaf.second.pop_back(); // 去除重复父结点
				if (root == (visitNum - 1)) { break; }
				vector<ID> lastChange(std::move(leaf.second));	// 备份到达该节点的访问路径
				int level = 0, vid = root > 0 ? (root + 1) : nodeNum; // 下次迭代从父结点的子节点开始（初始跳过第一个周期）
				stack<int> trace;
				trace.push(root);

				while (!trace.empty()) {	// 扩展该节点
					if (0 == (vid % nodeNum)) { ++vid; continue; }	//仓库不反转
					int curRoot = trace.top();
					ID p = curRoot / nodeNum, n = curRoot % nodeNum;
					if (curRoot == (visitNum - 1)) {	// 如果访问到末尾，弹出当前节点和其父节点
						trace.pop(); --level;
						visits[p][n] = 1 - visits[p][n];
						curRoot = trace.top();
						trace.pop(); --level;
						if (curRoot > root) {
							visits[curRoot / nodeNum][curRoot % nodeNum] = 1 - visits[curRoot / nodeNum][curRoot % nodeNum];
							vid = curRoot + 1;
						}
						else { break; }
					}
					else if (level == depth) {	// 达到深度，弹出当前节点
						visits[p][n] = 1 - visits[p][n];
						trace.pop(); --level;
						vid = curRoot + 1;
					}
					if (vid > curRoot) {		// 偏序搜索
						trace.push(vid);
						ID pid = vid / nodeNum, nid = vid % nodeNum;
						visits[pid][nid] = 1 - visits[pid][nid];
						if (visits[pid][nid]) { visits[pid][0] = 1; } // 客户变为1，仓库置1
						Price invCost = callModel(visits);		// 求改变访问后的库存
						if (invCost >= 0) {		//库存模型有解
							Price curTourCost = callLKH(sln, visits, pid);	// 求路由
							Price totalCost = sln.allTourCost + curTourCost - toursCost[pid] + invCost;
							cout << "cost after change : " << totalCost << ", best cost now : " << bestCost << endl;
							if (totalCost < bestCost) {	// 记录全局最优
								bestCost = totalCost;
								bestChange = lastChange;
								for (const auto &vid : trace._Get_container()) {
									bestChange.push_back(vid);
								}
							}	// 记录可扩展的叶子，不一定要topk
							if ((level + 1) == depth && vid < (visitNum - 1)) {
								leaf.first = totalCost;
								leaf.second = lastChange;
								for (const auto &v : trace._Get_container()) {
									leaf.second.push_back(v);
								}
								leaf.second.push_back(vid);	// 起始节点
								leaf.second.push_back(step + 1);// 为该叶子标注扩展周期
								if (topkLeaves.size() < leafNum) {
									topkLeaves.push_back(std::move(leaf));
								}
								else {
									sort(topkLeaves.begin(), topkLeaves.end());
									if (totalCost < topkLeaves.rbegin()->first) {
										topkLeaves[leafNum - 1] = std::move(leaf);
									}
								}
							}
						}
						++level; ++vid;
					}
				}
			}
			for (auto &l : topkLeaves) {
				leaves.push(std::move(l));
			}
		}
		sln.totalCost = bestCost;
		cout << "best result after iteration : " << bestCost << endl;
	}

	bool Solver::optimize(Solution &sln, ID workerId) {
		Log(LogSwitch::Szx::Framework) << "worker " << workerId << " starts." << endl;

		sln.init(input.periodnum(), input.vehicles_size(), Problem::MaxCost);
		iteratedModel(sln);
		treeSearch(sln, aux.m, aux.iterTime);

		Log(LogSwitch::Szx::Framework) << "worker " << workerId << " ends." << endl;
		return true;
	}

	void Solver::iteratedModel(Solution &sln) {
		ID nodeNum = input.nodes_size();
		ID vehicleNum = input.vehicles_size();
		const auto &vehicles(*input.mutable_vehicles());
		const auto &nodes(*input.mutable_nodes());

		MpSolver::Configuration mpCfg(MpSolver::InternalSolver::GurobiMip, env.timeoutInSecond(), true, false);
		MpSolver mp(mpCfg);

		// delivery[p, v, n] is the quantity delivered to node n at period p by vehicle v.
		Arr2D<Arr<Dvar>> delivery(input.periodnum(), vehicleNum, Arr<Dvar>(nodeNum));
		// x[p, v, n, m] is true if the edge from node n to node m is visited at period p by vehicle v.
		Arr2D<Arr2D<Dvar>> x(input.periodnum(), vehicleNum, Arr2D<Dvar>(nodeNum, nodeNum));

		// quantityLevel[n, p] is the rest quantity of node n at period p after the delivery and consumption have happened.
		Arr2D<Expr> quantityLevel(nodeNum, input.periodnum());
		// degrees[n] is the sum of in-coming edges.
		Arr2D<Arr<Expr>> degrees(input.periodnum(), vehicleNum, Arr<Expr>(nodeNum));

		// add decision variables.
		for (ID p = 0; p < input.periodnum(); ++p) {
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
			for (ID p = 0; p < input.periodnum(); ++p) {
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

		for (ID p = 0; p < input.periodnum(); ++p) {
			for (ID v = 0; v < vehicleNum; ++v) {
				Expr quantity;
				for (ID n = 0; n < nodeNum; ++n) {
					quantity += delivery[p][v][n];
				}
				// quantity matching constraint.
				mp.addConstraint(quantity == 0);
			}
		}

		for (ID p = 0; p < input.periodnum(); ++p) {
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
			for (ID p = 0; p < input.periodnum(); ++p) {
				holdingCost += (node.holdingcost() * quantityLevel.at(n, p));
			}
		}
		Expr routingCost;
		for (ID p = 0; p < input.periodnum(); ++p) {
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
			for (ID p = 0; p < input.periodnum(); ++p) {
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

		mp.setMipSlnEvent(subTourHandler);

		static const String TspCacheDir("TspCache/");
		System::makeSureDirExist(TspCacheDir);
		CachedTspSolver tspSolver(nodeNum, TspCacheDir + env.friendlyInstName() + ".csv");

		Solution curSln; // current solution.
		curSln.init(input.periodnum(), vehicleNum);
		auto nodeSetHandler = [&](MpSolver::MpEvent &e) {
			//if (e.getObj() > sln.totalCost) { e.stop(); } // there could be bad heuristic solutions.

			// OPTIMIZE[szx][0]: check the bound and only apply this to the optimal sln.

			lkh::CoordList2D coords; // OPTIMIZE[szx][3]: use adjacency matrix to avoid re-calculation and different rounding?
			coords.reserve(nodeNum);
			List<ID> nodeIdMap(nodeNum);
			List<bool> containNode(nodeNum);
			lkh::Tour tour;
			Expr nodeDiff;

			curSln.totalCost = 0, curSln.allTourCost = 0;
			for (ID p = 0; p < input.periodnum(); ++p) {
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
						curSln.allTourCost += aux.routingCost.at(*n, *m);
					}
				}
			}
			e.addLazy(nodeDiff >= 1);

			curSln.totalCost += e.getValue(holdingCost);
			if (curSln.totalCost < sln.totalCost) {
				Log(LogSwitch::Szx::Model) << "opt=" << curSln.totalCost << endl;
				cout << "allTourCost = " << curSln.allTourCost << endl;
				swap(curSln, sln);
			}
		};

		mp.setMipSlnEvent(nodeSetHandler);

		// solve.
		if (mp.optimize()) {
			// retrive solution.
			sln.totalCost = obj.getValue();

			for (ID p = 0; p < input.periodnum(); ++p) {
				auto &periodRoute(*sln.mutable_periodroutes(p));
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
		}
	}
#pragma endregion Solver

}
