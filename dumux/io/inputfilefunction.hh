#ifndef INPUT_FILE_FUNCTION_HH
#define INPUT_FILE_FUNCTION_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>

#include "../../dumux/external/csv.h"

// I replaced assert with if, since I don't know how to activate them (somewhere NDEBUG must be set)

namespace Dumux {

/**
 * Auxiliary class to pass a value or function in one argument via the input file, or via the grid file.
 * The type is chosen automatically.
 *
 * The type is determined from y, and x, that are read from the input file:
 *
 * Constant     returns one constant value
 *              y is set to a single value, x is not set
 *
 * Table        the value is calculated by linear interpolation from a table (y,x) given in the input file
 *              y and x are both set to multiple sampling points
 *              y and x are not set, but a parameter CSVFile = example.csv, is set in the same parameter group.
 *
 * Data         the value is taken from the grid file
 *              y and x are not set
 *
 * PerType      the value is constant and given per type (e.g. soil layer, root order), where type is taken from the grid file.
 *              y is set to multiple values, x is not set.

 * PerTypeIFF   the value is constant and given per type (e.g. soil layer, root order), where type is taken from another InputFileFunction calss
 *              y is set to multiple values, x is not set.
 *
 * TablePerType Multiple tables, named x0, x1, ... and y0, y1, ... The table number is chosen per type. (e.g. conductivities per root type)
 *              y and x are not set, x0,x1, ... , y0,y1, ... are set.
 *
 * All function types can be made sinusoidal over one day, by having the parameter Sinusoidal = True in the same parameter group.
 *
 * If grid data is used (data, perType, tablePerType), the data must be set with InputFileFunction::setGridData,
 * before calling InputFileFunction::f. Other data sources can be used using InputFileFunction::setGridData.
 */
class InputFileFunction {
public:

    enum {
        constant = 0, table = 1, data = 2, perType = 3, perTypeIFF = 4, tablePerType = 5, empty = -1
    };

    InputFileFunction() {  } // empty

    /**
     * Constructor 1: an input file function of all types
     *
     * call InputFileFunction::f(double x, size_t eIdx)
     *
     * @param groupName         name of the group in the .input file
     * @param nameY             name of the parameter
     * @param nameX             name of the holding sampling point coordinates (for table, tablePerType) (optional)
     * @param dataIdx           index, where the parameter values are located in the mesh file (for data) (optional)
     * @param typeIdx           index, where the type/order is located in the mesh file (for perType, talbePerType) (optional)
     * @þaram typeF             an input file function determining the type (for perType, talbePerType) (optional)
     */
    InputFileFunction(std::string groupName, std::string nameY, std::string nameX, int dataIdx, int typeIdx,
        InputFileFunction* typeF = nullptr) {

        nameY = groupName +"." + nameY; // full names
        nameX = groupName +"." + nameX; // full names
        dataIdx_ = dataIdx;
        typeIdx_ = typeIdx;
        nameY_ = nameY;
        sinusoidal_ = (hasParam(groupName +".Sinusoidal")) ? Dumux::getParam<bool>(groupName +".Sinusoidal") : false;
        if (Dumux::hasParam(nameY)) {
            yy_ = Dumux::getParam<std::vector<double>>(nameY);
            if (yy_.size() == 1) {
                type_ = constant;
            } else {
                if (Dumux::hasParam(nameX)) {
                    type_ = table;
                    xx_ = Dumux::getParam<std::vector<double>>(nameX);
                    table_.push_back( { xx_, yy_ });
                } else {
                    if (typeF != nullptr) {
                        type_ = perTypeIFF;
                        iff_ = typeF;
                    } else {
                        type_ = perType;
                    }
                }
            }
        } else {
            if (Dumux::hasParam(nameY + std::to_string(0))) {
                int m = 0;
                while (Dumux::hasParam(nameY + std::to_string(m))) { // count how many there are
                    m++;
                }
                for (int i = 0; i < m; i++) { // read them
                    xx_ = Dumux::getParam<std::vector<double>>(nameX + std::to_string(i));
                    yy_ = Dumux::getParam<std::vector<double>>(nameY + std::to_string(i));
                    table_.push_back( { xx_, yy_ });
                }
                type_ = tablePerType;
            } else {
                if  (hasParam(groupName +".CSVFile")) {
                    std::string name = Dumux::getParam<std::string>(groupName +".CSVFile");
                    std::string filename = getParam<std::string>(groupName +".CSVFile");
                    io::CSVReader<2> csv(filename);
                    csv.read_header(io::ignore_missing_column+io::ignore_extra_column, "x", "y");
                    std::vector<double> x, y;
                    double a,b;
                    while(csv.read_row(a,b)){
                        x.push_back(a);
                        y.push_back(b);
                    }
                    table_.push_back( { x, y });
                    type_ = table;
                } else {
                    type_ = data;
                }
            }
        }
        //cout();
    }

    /**
     * Constructor 2 (no tables):   an input file functions of types: constant, perType, perTypeIFF or data
     *                              i.e. no dependency of a variable
     *
     * a call to InputFileFunction::f(size_t eIdx) is sufficient
     */
    InputFileFunction(std::string groupName, std::string nameY, int dataIdx = 0, int typeIdx = 0, InputFileFunction* typeF = nullptr)
    :InputFileFunction(groupName, nameY, "no_valid_name", dataIdx, typeIdx, typeF) {

        //        assert(((type_==constant) || (type_==perType) || (type_==perTypeIFF) ||(type_==data)) &&
        //            "InputFileFunction: wrong type in constructor 2");
        if (!((type_==constant) || (type_==perType) || (type_==perTypeIFF) ||(type_==data))) {
            throw Dumux::ParameterException("InputFileFunction: wrong type in constructor 2, parsing "+
                nameY+", type = " + std::to_string(type_));
        }
    }

    /**
     * Constructor 3:   an input file functions of types constant, table, or perTypeIFF
     *                  i.e. all types without grid data
     *
     * a call to InputFileFunction::f(double x) is sufficient
     */
    InputFileFunction(std::string groupName, std::string nameY, std::string nameX, double defaultValue, InputFileFunction* typeF = nullptr)
    :InputFileFunction(groupName, nameY, nameX, 0, 0, typeF) {
        if (type_==data) {
            type_ = constant;
            yy_ = { defaultValue };
            std::cout << "\e[A"; // deletes last line (not in eclipse)
            cout();
        }
        //        assert(((type_==constant) || (type_==table) || (type_==perTypeIFF))  &&
        //            "InputFileFunction: wrong type in constructor 3");
        if (!((type_==constant) || (type_==table)|| (type_==perTypeIFF))) {
            throw Dumux::ParameterException("InputFileFunction: wrong type in constructor 3 parsing "+nameY+
                ", type = " + std::to_string(type_));
        }
    }

    /*
     * Constructor 4:    a hard coded a linear look up table
     *
     * it should not be necessary to use this constructor, a call to InputFileFunction::f(double x) is sufficient
     */
    InputFileFunction(std::vector<double> x, std::vector<double> y) {
        type_ = table;
        table_.push_back( { x, y });
        std::cout << "InputFileFunction: hard coded table \n";
    }

    /**
     * Scales the dependent variable with @param s in InputFileFunction::f calls (for unit conversions)
     */
    void setVariableScale(double s) {
        vs_ = s;
    }

    /**
     * Scales the function's return value with @param s in InputFileFunction::f calls (for unit conversions)
     */
    void setFunctionScale(double s) {
        fs_ = s;
    }

    /**
     * Evaluates the function
     *
     * @param x         the variable (optionally is scaled, @see setVariableScale)
     * @param eIdx      Dumux element index for grid data look up
     *
     * @return          the function return values (optionally is scaled, @see setFunctionScale)
     */
    double f(double x, size_t eIdx) const {

        double fs  = fs_;
        if (sinusoidal_)  {
            fs = fs * (sin((x/86400.)*2.* M_PI - 0.5*M_PI) + 1);
            // x is assumed to be in seconds, scaling later, e.g. for look up tables based on days
        }
        x = x*vs_; // scale variable

        switch (type_) {
        case constant: {
            return fs*yy_[0];
        }
        case table: {
			if(nameY_ != "Soil.IC.P")// to avoid error in case of small extrapolation
			{
				auto table = &table_[0];
				auto ip = x;
				const auto& range = table->first;
				const auto& values = table->second;

				// check bounds
				if (ip > range.back()) return values.back();
				if (ip < range[0]) return values[0];

				// if we are within bounds find the index of the lower bound
				const auto lookUpIndex = std::distance(range.begin(), std::lower_bound(range.begin(), range.end(), ip));
				if (lookUpIndex == 0)
					return values[0];
			}
            return fs*Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_[0]);
        }
        case data: {
            return fs*data_.at(eIdx);
        }
        case perType: {
            int t = int(data_.at(eIdx));
            if (t < 0) {
                std::cout << "InputFileFunction::perType: warning type is negative for element index " << eIdx <<": " << t << ", resuming with root order 0\n" << std::flush;
                t = 0;
            }
            return fs*yy_.at(t);
        }
        case perTypeIFF: {
            int t = int(iff_->f(x, eIdx));
            if (t < 0) {
                std::cout << "InputFileFunction::perTypeIFF: warning type is negative for element index " << eIdx <<": " << t << ", resuming with root order 0\n" << std::flush;
                t = 0;
            }
            return fs*yy_.at(t);
        }
        case tablePerType: {
            int t = int(data_.at(eIdx));
            if (t>=int(table_.size())) {
                std::cout << "InputFileFunction::tablePerType: warning InputFileFunction::f("+std::to_string(x)+", "
                    +std::to_string(eIdx)+"): "+nameY_+" table type > available tables, "+ std::to_string(t)+">="+std::to_string(table_.size()) +"\n";
                t = table_.size()-1;
            }
            if (t < 0) {
                std::cout << "InputFileFunction::tablePerType: warning type is negative for element index " << eIdx <<": " << t << ", resuming with root order 0\n" << std::flush;
                t = 0;
            }
            //            assert( t<table_.size() && "InputFileFunction::f: table type > available tables" );

            return fs*Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_.at(t));
        }
        default:
            throw Dumux::ParameterException("InputFileFunction: unknown function type");
        }
    }

    /**
     *
     */
    double f(size_t eIdx) const {
        //        assert((type_ != table &&  "InputFileFunction: call f(x) for table"));
        if ((type_ == table)||(type_ == tablePerType)) {
            throw Dumux::ParameterException("InputFileFunction::f(eIdx): call f(x) for function types with a variable (tables)");
        }
        return f(0, eIdx);
    }

    /**
     *
     */
    double f(double x) const {
        //        assert((type_ != data && type_ != perType && "InputFileFunction: call f(eIdx) for data or perType"));
        if ((type_ == data)||(type_ == perType)||(type_ == tablePerType)) {
            throw Dumux::ParameterException("InputFileFunction::f(x): call f(eIdx) for function types that depend on data");
        }
        return f(x, 0);
    }

    /**
     * Type index of the InputFileFunction
     */
    int type() const {
        return type_;
    }

    /**
     * Sets either parameter values per element index (for data), or
     * type values per element index (for perType, tablePerType)
     *
     * must be called before f(..) for data, perType, and tablePerType
     */
    void setData(const std::vector<double> d) {
        if ((type_ == data) || (type_ == perType) || (type_ == tablePerType)) { // otherwise, don't bother
            // std::cout << "InputFileFunction::setData for " << nameY_ << " \n";
            data_ = d;
        }
    }

    /**
     * Copies the data from the grid file using the values as
     * parameter values per element index (for data), or
     * type values per element index (for perType, tablePerType)
     *
     * must be called before f(..) for data, perType, and tablePerType
     */
    template<class GridData, class FVGridGeometry>
    void setGridData(const GridData& gridData, const FVGridGeometry& fvGridGeometry) {
        if ((type_ == data) || (type_ == perType) || (type_ == tablePerType)) { // otherwise, don't bother
            // std::cout << "InputFileFunction::setGridData for " << nameY_ << "\n";
            const auto& gridView = fvGridGeometry.gridView();
            data_.resize(gridView.size(0));
            const auto& elementMapper = fvGridGeometry.elementMapper();
            for (const auto& element : elements(gridView)) {
                const auto eIdx = elementMapper.index(element);
                const auto& elemParams = gridData.parameters(element);
                if (type_ == data) {
                    data_[eIdx] = elemParams[dataIdx_];
                } else { // type_ == perType
                    data_[eIdx] = elemParams[typeIdx_];
                }
            }
        }
    }

    /**
     * Quick info about the type of function
     */
    void cout() {
        std::string s = "";
        if (sinusoidal_) {
            s = ", sinusoidal over one day ";
        }
        switch (type_) {
        case(constant): std::cout << "InputFileFunction: "<< nameY_ << ": Constant"<< s << "\n"; break;
        case(table): std::cout << "InputFileFunction: "<< nameY_ << ": Table" << s << "\n";  break;
        case(data): std::cout << "InputFileFunction: "<< nameY_ << ": Data from Grid at index "<< dataIdx_ << s << "\n"; break;
        case(perType): std::cout << "InputFileFunction: "<< nameY_ << ": Constant per Type from Grid at index "<< typeIdx_ << s << "\n"; break;
        case(perTypeIFF): std::cout << "InputFileFunction: "<< nameY_ << ": Constant per Type from InputFileFunction" << s << "\n"; break;
        case(tablePerType): std::cout << "InputFileFunction: "<< nameY_ << ": Table per Type from Grid, " << table_.size()
                    << " types" << s << "\n"; break;
        default:
            std::cout << "InputFileFunction: unknown function type";
        }
    }

private:

    std::string nameY_;
    int type_ = empty;
    size_t dataIdx_ = -1;
    size_t typeIdx_ = -1;

    bool sinusoidal_=false;
    double vs_ = 1.; // for unit conversions
    double fs_ = 1.; // for unit conversions
    InputFileFunction* iff_ = nullptr;
    std::vector<double> xx_;
    std::vector<double> yy_;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> table_;
    std::vector<double> data_;

};

} // namespace DUMUX

#endif
