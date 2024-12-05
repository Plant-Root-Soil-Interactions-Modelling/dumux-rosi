from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
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
tid = 0

# Master node
if rank == 0:
    trials = Trials()
    eval_count = 0
    active_workers = size - 1  # All workers except the master

    while eval_count < max_evals or active_workers > 0:
        # Check for messages from workers
        status = MPI.Status()
        message = comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = status)

        if status.tag == 1:  # Worker is requesting a task
            if eval_count < max_evals:
                # Generate a new hyperparameter set
                param = fmin(
                    fn = lambda p: 0,  # Dummy function (we handle evaluations manually)
                    space = space,
                    algo = tpe.suggest,
                    max_evals = 1,
                    trials = trials,
                )
                # trials["tid"] = tid
                #

                comm.send(param, dest = status.source, tag = 2)
                eval_count += 1
            else:
                # No more tasks; tell the worker to stop
                comm.send(None, dest = status.source, tag = 0)

        elif status.tag == 3:  # Worker sending back results
            worker_result = pickle.loads(message)
            tid += 1
            doc = {
                "tid":tid,
                "result": worker_result["loss"],
                "spec": "tpe.suggest, max_evals 1",
                "misc": {"tid":tid },
                "state": worker_result["status"],
                "owner": "root",  # chatGPT told me
                "book_time": 0.,  # no glue
                "refresh_time": 0.,  # no glue
                "exp_key": "myexperiment",
                "cmd": None
                }

            print(worker_result)
            trials.insert_trial_doc(doc)
            trials.refresh()

            # Reduce active worker count if stopping
            if worker_result.get("stop", False):
                active_workers -= 1

    print("Optimization completed.")
    print("Best parameters:", trials.best_trial["result"])

# Worker nodes
else:
    while True:
        # Request a task from the master
        comm.send(None, dest = 0, tag = 1)  # Request a task
        task = comm.recv(source = 0, tag = MPI.ANY_TAG)

        if task is None:
            # No more tasks; terminate
            break

        # Evaluate the task
        result = objective(task)

        # Send the result back to the master
        comm.send(pickle.dumps(result), dest = 0, tag = 3)
