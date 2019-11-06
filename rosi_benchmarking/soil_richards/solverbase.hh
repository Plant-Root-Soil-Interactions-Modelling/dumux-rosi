



class SolverBase {

public:

    virtual void createGrid(boundsMin, boundsMax, numberOfCells);
    virtual void createGrid(std::string file);

    virtual void initialConditions()
    virtual void boundaryConditions()

    virtual void initialize(); //
    virtual void simulate(); //

    int pickIndex(pos); // vertex or element index
    double solutionAt(pos);

    virtual  getPoints(); // nodes
    virtual  getCells();


    double simTime = 0;
    std::vector<double> initialValues;
    std::vector<double> solution;

private:

};






void SolverBase::initialize() {

    // is the grid set (?)

}



