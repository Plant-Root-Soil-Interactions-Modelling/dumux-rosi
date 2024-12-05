from hyperopt import hp, Trials, STATUS_OK
from hyperopt.tpe import suggest
from mpi4py import MPI
import pickle

# MPI initialization
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# Define the objective function
def objective(params):
    x = params["x"]
    y = params["y"]
    loss = (x - 3) ** 2 + (y + 5) ** 2  # Quadratic loss function
    return {"loss": loss, "status": STATUS_OK}


# Define the search space
space = {
    "x": hp.uniform("x", -10, 10),
    "y": hp.uniform("y", -10, 10),
}

# Number of evaluations
max_evals = 50

# Master node
if rank == 0:
    trials = Trials()
    eval_count = 0
    active_workers = size - 1  # All workers except the master
    task_id = 0  # Unique task identifier

    while eval_count < max_evals or active_workers > 0:
        # Check for messages from workers
        status = MPI.Status()
        message = comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = status)

        if status.tag == 1:  # Worker is requesting a task
            if eval_count < max_evals:
                # Generate a new hyperparameter set with TPE
                new_trials = suggest(
                    n_startup_jobs = 10,  # Use random sampling for the first few trials
                    domain = None,
                    trials = trials,
                    max_evals = eval_count + 1,
                )
                trial = new_trials[-1]  # Get the new trial

                # Assign a unique task ID (tid)
                trial["tid"] = task_id
                task_id += 1

                # Send the task to the worker
                comm.send(pickle.dumps(trial), dest = status.source, tag = 2)
                eval_count += 1
            else:
                # No more tasks; tell the worker to stop
                comm.send(None, dest = status.source, tag = 0)

        elif status.tag == 3:  # Worker sending back results
            result = pickle.loads(message)
            trials.insert_trial_doc(result)
            trials.refresh()

            # Reduce active worker count if stopping
            if result.get("stop", False):
                active_workers -= 1

    print("Optimization completed.")
    print("Best parameters:", trials.best_trial["result"])

# Worker nodes
else:
    while True:
        # Request a task from the master
        comm.send(None, dest = 0, tag = 1)  # Request a task
        task_message = comm.recv(source = 0, tag = MPI.ANY_TAG)

        if task_message is None:
            # No more tasks; terminate
            break

        # Evaluate the task
        trial = pickle.loads(task_message)
        trial["result"] = objective(trial["misc"]["vals"])  # Compute the loss
        trial["state"] = STATUS_OK

        # Send the result back to the master
        comm.send(pickle.dumps(trial), dest = 0, tag = 3)
