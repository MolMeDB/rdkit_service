
#### Get running jobs
def run(sshClient, username, include_finished = False):
    if include_finished:
        jobs = sshClient.shell_exec("qstat -xwfu " + username, clear=False, sleepTime=5)
    else:
        jobs = sshClient.shell_exec("qstat -pwfu " + username, clear=False, sleepTime=5)

    jobs = parseRunningJobs([str(row).split() for row in jobs])

    return {"status": "ok", "total": len(jobs), "jobs": jobs}

def parseRunningJobs(jobs):
    result = []
    
    curr_job = False
    for row in jobs:
        if len(row) == 0:
            continue
        if str(row[0]).strip() == "Job":
            if curr_job:
                result.append(curr_job)
            curr_job = False
            if len(row) > 2 and str(row[1]) == "Id:":
                curr_job = {
                    "id": " ".join(row[2:]),
                }
            continue
        if curr_job != False and len(row) > 2 and str(row[1]).startswith("="):
            # Filter only MolMeDB jobs
            if row[0] == "Job_Name" and not str(row[2]).startswith("MMDB"):
                curr_job = False
                continue
            curr_job[row[0]] = " ".join(row[2:])
            
    if curr_job != False:
        result.append(curr_job)

    return result