//    std::cout << "Try to simulate a new Crootbox root system" << "\n" << std::flush;
//    auto rootSystem = std::make_shared<CRootBox::RootSystem>();
//    rootSystem->openFile(fileName);
//    rootSystem->initialize();
//    rootSystem->simulate(getParam<double>("RootSystem.Grid.DtInitial"));
//    rootSystem->write("rb_rootsystem.vtp");
//    grid = RootSystemGridFactory::makeGrid(*rootSystem);
//    std::cout << "created the thing \n" << "\n" << std::flush;

using Grid = GetPropType<TypeTag, Properties::Grid>;
// Dune::FoamGrid<1,3>;
// std::shared_ptr<Grid> grid;
auto fileName = getParam < std::string > ("RootSystem.Grid.File");

std::string dgf = ".dgf";
//if (std::equal(dgf.rbegin(), dgf.rend(), fileName.rbegin())) { // dgf
std::cout << "try to open dgf" << "\n" << std::flush;

auto periodic = std::bitset < 3 > ("110");
using GridView = GetPropType<TypeTag, Properties::GridView>;
using GlobalCoordinate = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
auto ll = GlobalCoordinate();
auto ur = GlobalCoordinate();
for (size_t i = 0; i < 3; i++) {
    ll[i] = -100;
    ur[i] = 100;
}
PeriodicNetworkGridManager<3> gridManager(ll, ur, periodic);
// GridManager<Grid> gridManager;

gridManager.init("RootSystem");// pass parameter group (see input file)
auto& grid = gridManager.grid();
const auto gridData = gridManager.getGridData();
//Grid* pgrid = &grid_;
//grid.reset(pgrid);
std::cout << "opened dgf\n" << "\n" << std::flush;
//}
