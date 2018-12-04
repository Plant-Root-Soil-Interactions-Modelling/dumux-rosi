#ifndef INPUT_FILE_FUNCTION_HH
#define INPUT_FILE_FUNCTION_HH

#include <dumux/common/math.hh>

namespace Dumux {

/**
 * Auxiliary class to pass a function via the input file, or via the grid file. The type is chosen automatically
 *
 * Constant:    returns one constant value
 *              y is set to one value
 *
 * Table:       the value is calculated by linear interpolaton
 *              y and x are set to the sampling points
 *
 * Data:        the value is taken from the grid file
 *              y and x are not set
 *
 * Per Type     values is given per Type (e.g. soil layer, root order), where type is taken from the grid type.
 *              y is set to multiple value (>= max(types)).
 *
 * Periodic     (optional). Periodic over one day, e.g. for root collar.
 *              y is set to exactly two values determining min and max
 *
 * If grid data is used (Data, or Per Type), the data must be set with InputFileFunction::setData.
 */
class InputFileFunction {
public:

    enum {
        constant = 0, table = 1, data = 2, perType = 3, periodic = 4
    };

    InputFileFunction(std::string nameY, std::string nameX, bool enablePeriodic = false) {

        try {
            yy_ = Dumux::getParam<std::vector<double>>(nameY);
            if (yy_.size() == 1) {
                type_ = constant;
                std::cout << "InputFileFunction: Constant " << nameY << "\n";
            } else {
                try {
                    xx_ = Dumux::getParam<std::vector<double>>(nameX);
                    type_ = table;
                    table_ = {xx_, yy_};
                    std::cout << "InputFileFunction: Table " << nameY << "\n";
                } catch (...) {
                    if ((enablePeriodic) && yy_.size() == 2) {
                        type_ = periodic;
                        std::cout << "InputFileFunction: Periodic " << nameY << "\n";
                    } else {
                        type_ = perType;
                        std::cout << "InputFileFunction: Constant per Type from Grid " << nameY << "\n";
                    }
                }
            }
        } catch (...) {
            type_ = data;
            std::cout << "InputFileFunction: Data from Grid " << nameY << "\n";
        }

    }

    double f(double x, size_t eIdx) {
        switch (type_) {
        case constant: {
            return yy_[0];
        }
        case table: {
            return Dumux::interpolate<Dumux::InterpolationPolicy::LinearTable>(x, table_);
        }
        case data: {
            return data_->at(eIdx);
        }
        case perType: {
            return yy_.at(size_t(data_->at(eIdx)));
        }
        case periodic: {
            return sin(x / (24. * 3600.) * 2. * M_PI) * (yy_[1] - yy_[0]) + 0.5 * (yy_[0] + yy_[1]);
        }
        default:
            throw Dumux::ParameterException("InputFileFunction: unknown function type");
        }
    }

    void setData(std::vector<double>* data) {
        data_ = data;
    }

private:
    int type_;
    std::vector<double> xx_;
    std::vector<double> yy_;
    std::pair<std::vector<double>, std::vector<double>> table_;
    std::vector<double>* data_ = nullptr;

};

} // namespace DUMUX

#endif
