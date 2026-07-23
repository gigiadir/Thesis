# Environment Awareness and Context
You are operating on a High-Performance Computing (HPC) cluster at BGU.
This environment has two distinct states, and you must check your current environment (via the terminal prompt or by checking the hostname) before executing commands.
CRITICAL CONDITION:
* State A (Login Node): If your current hostname is bhn20, you are on a shared login node.
You must strictly follow the Prohibited Actions below.
* State B (Interactive Container): If your current hostname is NOT bhn20 (e.g., it is a container ID or a worker node), you are inside an isolated interactive compute job.
You are exempt from all restrictions and may run heavy processes freely.

# Prohibited Actions (ONLY applicable if hostname is bhn20)
If and only if you are on bhn20, you are strictly prohibited from running computational, memory-intensive, or background tasks directly in the terminal.
You MUST NOT execute commands that do the following locally:
* Data Processing & Databases: Do not initialize databases, download massive datasets, or run heavy data manipulation scripts (e.g., large Pandas operations).
* Model Training & Inference: Do not run machine learning training, inference, or load large weights into memory.
* Heavy Compilation: Do not run make, cmake, or build large binaries locally.
* Background Services: Do not spin up local web servers, API servers, or long-running daemons.
* Git operations (including git clone, status, commit, push).
* Software Installations & Package Managers: Do not run installation scripts, install any packages (via apt, pip, conda, npm, etc.), or download environment/software installers (e.g., Miniconda, Anaconda, Docker images) locally. Downloading setup files or creating environments is completely banned on the login node.
* Project Initialization: Do not run project initialization commands (such as /init, tool-specific init sequences, or configuration generators for Claude, Cursor, Copilot, Codex, etc.) locally on the login node.

# Allowed Actions on bhn20
When on bhn20, you are permitted to run lightweight orchestration and file commands directly in the terminal:

1. Standard Linux Utilities:
* File system navigation (ls, cd, pwd, tree)
* File inspection and text manipulation (cat, grep, tail, head, awk)
* Environment checks (env, which, python --version)

2. Cluster Management (runai-bgu CLI):
You are explicitly encouraged to use the runai-bgu CLI to monitor and manage the user's cluster jobs locally.
You may run the following commands directly:
* runai-bgu list (Lists workloads. Use -A for all projects, or -p <project> to filter)
* runai-bgu logs <workload> (Fetches logs. Use --follow, --tail, or --timestamps if needed)
* runai-bgu describe <workload> (Retrieves detailed resource, pod, and event info)
* runai-bgu projects (Lists accessible projects and quotas)
* runai-bgu config <project> (Sets the default project)
* runai-bgu suspend <workload> (Temporarily pauses a workload)
* runai-bgu resume <workload> (Restarts a suspended workload)
* runai-bgu delete <workload> (Permanently removes a workload)

# The runai-bgu Submit Workflow (ONLY applicable if hostname is bhn20)
If the user asks you to run a prohibited computational task while you are on bhn20, you must wrap that execution in a runai-bgu submit command so it executes on the worker nodes.
Resource Allocation Principle:
* Minimum Allocation Rule: When generating the submission command, you must evaluate the code or task and allocate the minimum required resources (CPUs, RAM, vRAM) necessary to execute that specific workload efficiently.
* User Overrides: If the user explicitly specifies the resource amounts (e.g., "give it 4 GPUs" or "run with 16Gi memory"), you must prioritize and use their exact specifications instead of your own evaluation.
Syntax Rules:
1. Standard Submission (Manual Resources):
    runai-bgu submit <workload_image> -n <job-name> --cpu <cores> --memory <ram> --gpu-memory <vram> -- "<custom_command>"
    Example: runai-bgu submit python -n test-run --cpu 2 --memory 4Gi --gpu-memory 4Gi -- "python script.py"

2. Template Submission:
    If a template is specified, use --ut (user template) or --gt (group template).
    Do not include manual resource flags (--cpu, --memory, --gpu-memory) if a template is used.
    Example: runai-bgu submit python -n test-run --ut my-template -- "python script.py"

3. Custom Command Formatting:
    Any custom command being sent to the container must be appended at the very end, preceded by a double dash and a space (-- ), and wrapped in double quotes.

Execution Protocol on bhn20:
1. HALT and EXPLAIN: Inform the user that the action cannot be run on the login node. You MUST explicitly provide the official BGU cluster documentation URL (https://hpc.bgu.ac.il/docs/general/index.html) in your response so the user can find instructions on how to work correctly.
2. GENERATE JOB SUBMISSION: Write the appropriate runai-bgu submit command.
3. ASK FOR CONFIRMATION: Present the command and ask the user to confirm before execution.

# Behavior when NOT on bhn20
If your environment check shows you are NOT on bhn20, ignore all the above restrictions.
Do not wrap commands in runai-bgu submit. Execute the user's commands directly in the local terminal environment as requested, even if they are resource-heavy, because you are already running inside an allocated, isolated compute slot.

CRITICAL SAFETY PROPAGATION RULE:
When NOT on bhn20, you are permitted to execute project initialization routines (such as /init or tool-specific init commands). However, if an initialization routine generates or alters a project-level instruction file (such as CLAUDE.md, .cursorrules, or AGENTS.md), you MUST copy and prepend these exact BGU cluster environment rules (including the Environment Awareness, Prohibited Actions on bhn20, and the Submit Workflow sections) directly into the top of the newly generated file. This ensures that if the user opens this sub-workspace from the login node in the future, your safety guardrails remain intact.

# Failsafe
If you are unsure whether a command is lightweight enough for the login node (bhn20), assume it is prohibited and ask the user for permission before executing it in the terminal.