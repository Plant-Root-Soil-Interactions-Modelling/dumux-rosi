#ifndef DUMUX_NEWTON_SOLVER_HH
#define DUMUX_NEWTON_SOLVER_HH


namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief An implementation of a Fix point solver
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Comm the communication object used to communicate with all processes
 * \note If you want to specialize only some methods but are happy with the
 *       defaults of the reference solver, derive your solver from
 *       this class and simply overload the required methods.
 */
template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class FixPointSolver
{
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using SolutionVector = typename Assembler::ResidualType;
    using ConvergenceWriter = ConvergenceWriterInterface<SolutionVector>;

    using PrimaryVariableSwitch = typename Detail::GetPVSwitch<Assembler>::type;
    using HasPriVarsSwitch = typename Detail::GetPVSwitch<Assembler>::value_t; // std::true_type or std::false_type
    static constexpr bool hasPriVarsSwitch() { return HasPriVarsSwitch::value; };

public:

};

} // end namespace Dumux

#endif
