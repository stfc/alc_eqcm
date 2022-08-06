!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that automatically sets up the submission script 
! directives for High-Performance-Computing (HPC)
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
! 
! Author    - i.scivetti May 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module hpc_setup
  Use constants,        Only : 
  Use fileset,          Only : file_type, &
                               FILE_SET_EQCM, &
                               FILE_HPC_SETTINGS
  Use input_types,      Only : in_integer, &
                               in_integer_array,   & 
                               in_logic,   & 
                               in_string
  Use numprec,          Only : wi, &
                               wp
  Use process_data,     Only : capital_to_lower_case, &
                               get_word_length
  Use references,       Only : web_slurm
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  ! Module elements
  Type :: module_type
    Integer(Kind=wi)   :: num=0
    Character(Len=256) :: element(100)
  End Type

  ! Type for the modelling related variables 
  Type, Public :: hpc_type
    Private
    ! Flag to generate simulation files
    Logical,            Public :: generate=.False.
    ! Script file name
    Character(Len=256), Public :: script_name
    ! Number of processes 
    Type(in_integer),  Public :: number_mpi_tasks 
    ! Number of nodes 
    Type(in_integer),  Public :: number_nodes
    ! Number of CPUs per node
    Type(in_integer),  Public :: cpus_per_node 
    ! Threads per process
    Type(in_integer),  Public :: threads_per_process
    ! Processes/tasks per node
    Integer(Kind=wi)          :: processes_per_node
    ! Total CPUs 
    Integer(Kind=wi),  Public :: total_cpus
    ! HPC machine name
    Type(in_string),   Public :: machine_name
    ! project name
    Type(in_string),   Public :: project_name
    ! job name
    Type(in_string),   Public :: job_name
    ! Plataform
    Type(in_string),   Public :: distribution
    ! Plataform
    Type(in_string),   Public :: platform
    ! parallelism_type
    Type(in_string),   Public :: parallelism_type
    ! queue
    Type(in_string),   Public :: queue
    ! memory_per_cpu (in MB)
    Type(in_integer),  Public :: memory_per_cpu
    ! required modules
    Type(in_string),   Public :: modules_info 
    Type(module_type), Public :: modules
    ! Executable
    Type(in_string),   Public :: executable
    ! time_limit 
    Type(in_integer_array),  Public :: time_limit
    ! Submission command
    Type(in_string),   Public :: exec_options    
    ! MKL 
    Type(in_logic),    Public :: mkl    
  Contains
     Private
     Procedure, Public  :: init_input_variables   =>  allocate_input_variables
     Final              :: cleanup
  End Type hpc_type

  Public :: check_hpc_settings, summary_hpc_settings  
  Public :: build_hpc_script

Contains

  Subroutine allocate_input_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential variables to build scripts for simulation 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(hpc_type), Intent(InOut)  :: T

    Integer(Kind=wi)   :: fail(1)
    Character(Len=256) :: message

    Allocate(T%time_limit%value(3),           Stat=fail(1))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of HPC variables for buiding scripts for job submission'
      Call error_stop(message)
    End If

  End Subroutine allocate_input_variables


  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(hpc_type) :: T

    If (Allocated(T%time_limit%value)) Then
      Deallocate(T%time_limit%value)
    End If

  End Subroutine cleanup 

  Subroutine check_hpc_settings(files, code_format, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check input HPC setings via &block_HPC_settings
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Character(Len=*),  Intent(In   ) :: code_format
    Type(hpc_type),    Intent(InOut) :: hpc_data

    Character(Len=256) :: message
    Character(Len=256) :: messages(3)
    Character(Len=256) :: error_block, aux
    Character(Len=2)   :: par
    Integer(Kind=wi)   :: i, ic
    Logical :: error

    par = '()'

    error_block = '***ERROR in &Block_hpc_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    Call info(' ', 1)
    Write (messages(1),'(1x,a)')  '***IMPORTANT: by specification of "&block_hpc_settings", the user has'
    Write (messages(2),'(1x,a)')  'requested to generate script files for HPC submission.'
    Call info(messages, 2) 

    !  Machine
    If (.Not. hpc_data%machine_name%fread) Then
      Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify the name of the HPC machine via&
                                & "machine_name".'
      Call error_stop(message)
    End If

    ! Platform
    If (hpc_data%platform%fread) Then
      If (Trim(hpc_data%platform%type) /= 'slurm') Then
         Write (message,'(2(1x,a))') Trim(error_block), &
                            &'Invalid specification for directive "platform". To date, only  '
         Call error_stop(message)
      End If
    End If

    !  Job name
    If (.Not. hpc_data%job_name%fread) Then
      Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify the name of the HPC job via&
                                & "job_name".'
      Call error_stop(message)
    Else
      error=.False. 
      i=1
      ! Find if there is a parenthesis
      Do While (i <= Len_Trim(hpc_data%job_name%type) .And. (.Not. error))
        ic = Index(par, hpc_data%job_name%type(i:i))
        If (ic > 0) Then
          Write (message,'(2(1x,a))') Trim(error_block), 'Specification for "job_name" must NOT contain parenthesis.&
                                    & Please change.'    
          Call error_stop(message)
        End If
        i=i+1
      End Do
    End If

    If (hpc_data%project_name%fread) Then
      error=.False.
      i=1
      ! Find if there is a parenthesis
      Do While (i <= Len_Trim(hpc_data%project_name%type) .And. (.Not. error))
        ic = Index(par, hpc_data%project_name%type(i:i))
        If (ic > 0) Then
          Write (message,'(2(1x,a))') Trim(error_block), 'Specification for "project_name" must NOT contain parenthesis.&
                                    & Please change.'
          Call error_stop(message)
        End If
        i=i+1
      End Do      
    End If

    ! Total number of processes/tasks
    If (hpc_data%number_mpi_tasks%fread) Then
      If (hpc_data%number_mpi_tasks%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "number_mpi_tasks" directive'
        Call error_stop(message)
      Else
        If (hpc_data%number_mpi_tasks%value < 1) Then
          Write (message,'(1x,2a)') Trim(error_block), ' Value for "number_mpi_tasks" must be equal to 1 at least'
          Call error_stop(message)
        End If
      End If
    Else 
      Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify directive "number_mpi_tasks"'
      Call error_stop(message)
    End If

    ! Total number of nodes 
    If (hpc_data%number_nodes%fread) Then
      If (hpc_data%number_nodes%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "number_nodes" directive'
        Call error_stop(message)
      Else
        If (hpc_data%number_nodes%value < 1) Then
          Write (message,'(1x,2a)') Trim(error_block), ' Value for "number_nodes" must be equal to 1 at least'
          Call error_stop(message)
        End If
      End If
    End If

    ! CPUs per node 
    If (hpc_data%cpus_per_node%fread) Then
      If (hpc_data%cpus_per_node%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "cpus_per_nodes" directive'
        Call error_stop(message)
      Else
        If (hpc_data%cpus_per_node%value < 1) Then
          Write (message,'(1x,2a)') Trim(error_block), ' Value for "cpus_per_node" must be equal to 1, at least'
          Call error_stop(message)
        End If
      End If
    End If

    ! parallelism_type
    If (hpc_data%parallelism_type%fread) Then
      If (Trim(hpc_data%parallelism_type%type) /= 'mpi-only' .And. &
         Trim(hpc_data%parallelism_type%type) /= 'mpi-openmp') Then
        Write (message,'(2(1x,a))') Trim(error_block), &
                                  &'Invalid specification for directive "parallelism_type". Options are MPI-only and MPI-OpenMP'
        Call error_stop(message)
      End If
    Else
      Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify directive "parallelism_type".&
                                & Options are: MPI and MPI-OpenMP' 
      Call error_stop(message)
    End If

    ! Threads per process 
    If (hpc_data%threads_per_process%fread) Then
      If (hpc_data%threads_per_process%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "threads_per_processs" directive'
        Call error_stop(message)
      Else
        If (hpc_data%threads_per_process%value < 1) Then
          Write (message,'(1x,2a)') Trim(error_block), ' Value for "threads_per_process" must be equal to 1, at least'
          Call error_stop(message)
        End If
        If (Trim(hpc_data%parallelism_type%type) == 'mpi-only') Then
          Write (message,'(2(1x,a))') Trim(error_block), 'Specification of directive "threads_per_process" is incompatible&
                                  & with MPI-only. Please remove it'
          Call error_stop(message)
        End If
      End If
    Else 
      If (Trim(hpc_data%parallelism_type%type) == 'mpi-openmp') Then 
        Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify directive "threads_per_process"'
        Call error_stop(message)
      Else
        hpc_data%threads_per_process%value=1
      End If
    End If

    ! queue
    If (.Not. hpc_data%queue%fread) Then
      Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify the partition to run&
                                & the simulation via directive "queue"'
      Call error_stop(message)
    End If

    ! time limit for the simulation
    If (hpc_data%time_limit%fread) Then
      If (hpc_data%time_limit%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "time_limit" directive.&
                                 & Format must be "  days  hours  minutes  " (all interger numbers)'
        Call error_stop(message)
      Else
        Do i =1, 3
          If (hpc_data%time_limit%value(i) < 0) Then
            Write (message,'(1x,2a)') Trim(error_block), ' None of the values provided in "time_limit" must be&
                                    & lower than zero.'
            Call error_stop(message)
          End If
        End Do 
        If (hpc_data%time_limit%value(1) == 0 .And. &
           hpc_data%time_limit%value(2) == 0 .And. &
           hpc_data%time_limit%value(3) == 0 ) Then
           Write (message,'(1x,2a)') Trim(error_block), ' Limit for the simulation time (time_limit) has been set&
                                   & to zero! Please change'     
           Call error_stop(message)
        End If
        If (hpc_data%time_limit%value(3) > 59) Then
          Write (message,'(1x,2a)') Trim(error_block), ' The amount of minutes cannot exceed 59.&
                                   & Please change'
          
          Call error_stop(message)
        End If
        If (hpc_data%time_limit%value(3) > 0 .And. hpc_data%time_limit%value(2) >= 24) Then
          Write (messages(1),'(1x,2a)') Trim(error_block), ' Wrong format for "time_limit". A day cannot exceed 24 hours.'
          Write (messages(2),'(1x,a)')  'Maximum allowed time for a day is 23 hours and 59 minutes.'
          Write (messages(3),'(1x,a)')  'If more time is needed, please add days, bearing in mind the time constraints for&
                                       & the HPC available.' 
          Call info(messages, 3)
          Call error_stop(' ')
        End If
      End If
    Else
      Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify directive "time_limit"'
      Call error_stop(message)
    End If 

    ! Memory per CPU
    If (hpc_data%memory_per_cpu%fread) Then
      If (hpc_data%memory_per_cpu%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "memory_per_cpu"&
                                 & directive.'
        Call error_stop(message)
      Else
        If (hpc_data%memory_per_cpu%value <= 0) Then
          Write (message,'(1x,2a)') Trim(error_block), ' Value of "memory_per_cpu" must be positive.'
          Call error_stop(message)
        End If
      End If
    End If

    ! Executable options 
    If (hpc_data%exec_options%fread) Then
      If (hpc_data%exec_options%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) specification for "exec_options"&
                         & directive.'
        Call error_stop(message) 
      End If
      error=.False.
      ! Check is inappropriate specification has been set in "exec_options" 
      If (Index(hpc_data%exec_options%type, 'mpirun') > 0 ) Then
        aux='mpirun'
        error=.True. 
      Else If (Index(hpc_data%exec_options%type, ' -n ') > 0 ) Then
        aux='-n'
        error=.True. 
      Else If (Index(hpc_data%exec_options%type, ' -np ') > 0 ) Then
        aux='-np'
        error=.True. 
      Else If (Index(hpc_data%exec_options%type, 'srun') > 0 ) Then
        aux='srun'
        error=.True. 
      End If
      If (error) Then
        Write(message, '(4a)') Trim(error_block), ' Specification for directive "exec_options" must not include "',&
                             & Trim(aux), '". Please review this setting and rerun.' 
        Call error_stop(message)
      End If 
    End If

    If (hpc_data%mkl%fread) Then
      If (hpc_data%mkl%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) specification for "mkl"&
                 & directive.'
      End  If
    Else
      hpc_data%mkl%stat=.False. 
    End If

    ! Executable  
    If (hpc_data%executable%fread) Then
      If (hpc_data%executable%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) specification for "executable"&
                         & directive.'
        Call error_stop(message) 
      End If
    Else 
      Write (message,'(2(1x,a))') Trim(error_block), 'The user must specify the executable of the code via&
                                & directive "executable". This setting can also be specified as a path&
                                & (including the executable).'
      Call error_stop(message)
    End If

    ! Check if the code and the parallelism_type are compatible
    aux=code_format
    Call capital_to_lower_case(aux) 
    If (Trim(aux)=='vasp') Then
      If (Trim(hpc_data%parallelism_type%type) /= 'mpi-only') Then
        Write (message,'(2(1x,a))') Trim(error_block), 'VASP can only be set for execution using option "mpi-only" for&
                                  & directive "parallelism_type". Please change.'
        Call error_stop(message)
      End If
    End If

    If (Trim(hpc_data%parallelism_type%type) == 'mpi-only' .And. hpc_data%threads_per_process%value > 1) Then
      Write (message,'(2(1x,a))') Trim(error_block), 'MPI-only jobs are not compatible with values of "threads_per_process"&
                                & larger than 1. Please change.'
      Call error_stop(message)
    End If
   
    ! Check SCARF related settings
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (Trim(hpc_data%machine_name%type)=='scarf') Then
      ! Platform
      If (Trim(hpc_data%platform%type) /= 'slurm') Then
        hpc_data%platform%type='slurm'
        Write (message,'(1x,a)') '***IMPORTANT: SCARF only allows the SLURM platform. ALC_EQCM will&
                               & change the assigned directive'
      End If
    ! Queue (check if it is a valid partition) 
      If (Trim(hpc_data%queue%type) /= 'ibis' .And. &
         Trim(hpc_data%queue%type) /= 'magnacarta' .And. &
         Trim(hpc_data%queue%type) /= 'derevolutionibus' .And. &
         Trim(hpc_data%queue%type) /= 'demagnete' .And. &
         Trim(hpc_data%queue%type) /= 'numanlys-cpu' .And. &
         Trim(hpc_data%queue%type) /= 'fbioctopus' .And. &
         Trim(hpc_data%queue%type) /= 'devel' .And. &
         Trim(hpc_data%queue%type) /= 'scarf' .And. &
         Trim(hpc_data%queue%type) /= 'gpu' .And. &
         Trim(hpc_data%queue%type) /= 'gpu-exclusive' .And. &
         Trim(hpc_data%queue%type) /= 'preemptable') Then
         Write (message,'(1x,4a)') Trim(error_block), &
                           &'Provided option "', Trim(hpc_data%queue%type), '" for directive "queue"&
                           & (partition) does not exist in SCARF'
         Call error_stop(message)
      End If

      error=.False. 
      If (hpc_data%time_limit%value(1) == 7) Then
        If (hpc_data%time_limit%value(2) /= 0 .Or. hpc_data%time_limit%value(3) /= 0) Then
          error=.True.
        End If
      End If
      If (hpc_data%time_limit%value(1) > 7) Then
        error=.True.
      End If
 
      If (error) Then
        Write (message,'(1x,2a)') Trim(error_block), ' Time limit for jobs in SCARF cannot exceed 7 days. Please set a&
                           & lower amount using directive "time_limit"'
        Call error_stop(message)
      End If

    Else
      Write (message,'(1x,4a)') Trim(error_block), ' The user must specify the "platform" for the machine ',&
                             & Trim(hpc_data%machine_name%type), '. To date, SLURM is the only option.'
      Call error_stop(message)
    End If 

    ! Total number of CPUs
    hpc_data%total_cpus=hpc_data%number_mpi_tasks%value*hpc_data%threads_per_process%value

    If (hpc_data%cpus_per_node%fread) Then
      If ( .Not. hpc_data%number_nodes%fread) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Directive "CPUs_per_node" is insuficient without the specification&
                                  & of "number_nodes". Please review the settings and rerun.'
        Call error_stop(message)
      End If
    End If

    ! Determine the tasks per node
    If (hpc_data%number_nodes%fread .And. hpc_data%cpus_per_node%fread) Then
      If (mod(hpc_data%number_mpi_tasks%value,hpc_data%number_nodes%value) /= 0)then
        Write (message,'(2(1x,a))') Trim(error_block), 'The ratio between "number_mpi_tasks" and "number_nodes"&
                                  & leads to a non-zero remainder. Please change to optimise the distribution.'
        Call error_stop(message)
      Else
        hpc_data%processes_per_node=hpc_data%number_mpi_tasks%value/hpc_data%number_nodes%value
        If (hpc_data%processes_per_node*hpc_data%threads_per_process%value /= hpc_data%cpus_per_node%value) Then
          If (Trim(hpc_data%parallelism_type%type) == 'mpi-only') Then
            Write (messages(1),'(1x,a,2(a,i4),a)') Trim(error_block), ' The resulting number of operations per node is ', &
                                    & hpc_data%processes_per_node, ', which is different from the specified ',& 
                                    & hpc_data%cpus_per_node%value, ' CPUs per node.'
            Write (messages(2),'(1x,a)') 'This configuration does not optimise the computational resources. Please change the&
                                    & settings for "number_mpi_tasks", "number_nodes" and/or "cpus_per_node".'
  
            Call info(messages, 2)
            Call error_stop(' ')
          ElseIf (Trim(hpc_data%parallelism_type%type) == 'mpi-openmp') Then
            Write (messages(1),'(1x,a,2(a,i4),a)') Trim(error_block), ' The resulting number of operations&
                                    & per node is ', hpc_data%processes_per_node, ', which is different from the specified ',& 
                                    & hpc_data%cpus_per_node%value, ' CPUs per node.'
            Write (messages(2),'(1x,a)') 'This configuration does not optimise the computational resources. Please change the&
                                    & settings for "number_mpi_tasks", "number_nodes", "cpus_per_node" and/or&
                                    & "threads_per_process".'
  
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End If
      End If
    End If

  End Subroutine check_hpc_settings

  Subroutine build_hpc_script(files, code_format, hpc_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Build HPC scripts to submit simulations
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Character(Len=*),  Intent(In   ) :: code_format
    Type(hpc_type),    Intent(InOut) :: hpc_data

    Integer(Kind=wi)   :: iunit, i
    Character(Len=256) :: timeformat, aux, filein, fileout

    hpc_data%script_name='hpc_script-'//Trim(code_format)//'.sh'

    ! Open FILE_SET_SIMULATION file
    Open(Newunit=files(FILE_HPC_SETTINGS)%unit_no, File=files(FILE_HPC_SETTINGS)%filename,Status='Replace')
    iunit=files(FILE_HPC_SETTINGS)%unit_no

    ! Get time from input time_limit
    Call get_hpc_time(hpc_data%time_limit%value, timeformat)

    If (Trim(hpc_data%platform%type) == 'slurm') Then
      Write (iunit,'(a)')  '#!/bin/bash'
      Write (iunit,'(2a)') '## Submission script for ',  Trim(hpc_data%machine_name%type)
      If (hpc_data%project_name%fread) Then
        Write (iunit,'(3a)') '#SBATCH --comment=',   Trim(hpc_data%project_name%type),   '    # Project name'
      End If
      Write (iunit,'(3a)') '#SBATCH --job-name=',  Trim(hpc_data%job_name%type),         '    # job name'
      Write (iunit,'(a)') '#SBATCH -o %J.out'
      Write (iunit,'(a)') '#SBATCH -e %J.err'
      If (hpc_data%time_limit%value(1)>0) Then
        Write (iunit,'(3a)') '#SBATCH --time=',    Trim(timeformat),          '        # days-hh:mm:ss'
      Else
        Write (iunit,'(3a)') '#SBATCH --time=',    Trim(timeformat),          '        # hh:mm:ss'
      End If
      Write (iunit,'(a)')  '#'
      If (hpc_data%memory_per_cpu%fread) Then
        Write (aux, *) hpc_data%memory_per_cpu%value
        Write (iunit,'(3a)') '#SBATCH --mem-per-cpu=',  Trim(Adjustl(aux)) , '  # Megabytes'
      End If
      Write (iunit,'(3a)')     '#SBATCH --partition=', Trim(hpc_data%queue%type),      '    # queue (partition)' 
      Write (aux, *) hpc_data%number_mpi_tasks%value
      Write (iunit,'(2a)') '#SBATCH --ntasks=', Adjustl(Trim(aux))
      If (hpc_data%number_nodes%fread) Then
        Write (aux, *) hpc_data%number_nodes%value
        Write (iunit,'(2a)') '#SBATCH --nodes=', Adjustl(Trim(aux))
        If (hpc_data%cpus_per_node%fread) Then
          Write (aux, *) hpc_data%processes_per_node
          Write (iunit,'(2a)') '#SBATCH --ntasks-per-node=', Adjustl(Trim(aux))
        End If 
      End If 
      If (Trim(hpc_data%parallelism_type%type) == 'mpi-openmp') Then      
        Write (aux, *) hpc_data%threads_per_process%value
        Write (iunit,'(2a)') '#SBATCH --cpus-per-task=', Adjustl(Trim(aux))
        Write (iunit,'(a)') ' '
        Write (iunit,'(2a)')   'export OMP_NUM_THREADS=', Adjustl(Trim(aux))
      End If
      If (hpc_data%mkl%stat) Then
        If (Trim(hpc_data%parallelism_type%type) == 'mpi-openmp') Then      
          Write (iunit,'(2a)') 'export MKL_NUM_THREADS=', Adjustl(Trim(aux))
        Else
          Write (iunit,'(a)')  'export MKL_NUM_THREADS=1'
        End If
      End If

      Write (iunit,'(a)') ' '
      Write (iunit,'(a)') '## Load required modules'
      Do i=1, hpc_data%modules%num
        Write (iunit,'(2a)') 'module load ', Trim(hpc_data%modules%element(i))
      End Do
      Write (iunit,'(a)') ' '
      Write (iunit,'(a)') '## Define executable'
      Write (iunit,'(3a)') 'exec="', Trim(hpc_data%executable%type), '"'
      Write (iunit,'(a)') ' '
      If (hpc_data%exec_options%fread) Then 
        Write (aux,*) hpc_data%number_mpi_tasks%value
        Write (aux,*) Trim(aux)//' '//Trim(hpc_data%exec_options%type) 
      Else
        Write (aux,*) hpc_data%number_mpi_tasks%value
      End If
      Write (iunit,'(a)') '## Execute job'
      If (Trim(code_format)=='vasp') Then
        Write (iunit,'(3a)') 'mpirun -np ', Trim(Adjustl(aux)), ' $exec'
      Else If (Trim(code_format)=='cp2k') Then
        filein ='input.'//Trim(code_format)
        fileout='output.'//Trim(code_format)
        Write (iunit,'(6a)') 'mpirun -np ', Trim(Adjustl(aux)), ' $exec -i ', Trim(filein), ' -o ', Trim(fileout)
      Else If (Trim(code_format)=='castep') Then  
        Write (iunit,'(3a)') 'mpirun -np ', Trim(Adjustl(aux)), ' $exec model' 
      Else If (Trim(code_format)=='onetep') Then  
        Write (iunit,'(3a)') 'mpirun -np ', Trim(Adjustl(aux)), ' $exec  model.dat > model.out' 
      End If
    End If

    close(iunit)

 
  End Subroutine build_hpc_script

  Subroutine summary_hpc_settings(hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print summary of HPC setings to inform the user
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(hpc_type),    Intent(In   ) :: hpc_data

    Character(Len=256)  :: messages(16)
    Real(Kind=wp)       :: cpu_hours

    Call info(' =================================', 1)
    Call info(' Summary of HPC settings (per job)', 1)
    Call info(' =================================', 1)
    Write (messages(1),'(1x,2a)')   '- job will be executed using ', Trim(hpc_data%parallelism_type%type)
    Write (messages(2),'(1x,a,i5)') '- requested processes/tasks: ', hpc_data%number_mpi_tasks%value
    Call info(messages,2)

    If (hpc_data%number_nodes%fread) Then
      Write (messages(1),'(1x,a,i5)') '- requested nodes:           ', hpc_data%number_nodes%value 
      Call info(messages,1)
    End If

    If (hpc_data%cpus_per_node%fread) Then
      Write (messages(1),'(1x,a,i5)') '- CPUs per nodes:            ', hpc_data%cpus_per_node%value 
      Call info(messages,1)
    End If

    If (hpc_data%mkl%stat) Then
      Write (messages(1),'(1x,a)')    '- MKL libraries will be used for the job execution'
      Call info(messages,1)
    End If
   
    If (Trim(hpc_data%parallelism_type%type) == 'mpi-openmp') Then
      Write (messages(1),'(1x,a,i5)') '- threads per process/task:  ', hpc_data%threads_per_process%value
      Call info(messages,1)
    End If
    Write (messages(1),'(1x,a,i5)') '- total number of CPUs:      ', hpc_data%total_cpus
    cpu_hours=hpc_data%total_cpus*(hpc_data%time_limit%value(1)*24.0_wp + hpc_data%time_limit%value(2)*1.0_wp +&
                                & hpc_data%time_limit%value(3)/60.0_wp)
    Write (messages(2),'(1x,a,f9.2)') '- maximum CPU.Hours booked: ', cpu_hours
    If (hpc_data%memory_per_cpu%fread) Then
      Write (messages(3),'(1x,a,i7)') '- total memory (MBytes):   ', hpc_data%total_cpus*hpc_data%memory_per_cpu%value
    Else
      Write (messages(3),'(1x,a)')    '- the total memory will be determined by the allocated nodes'
    End If 
    Call info(messages,3)

    Call info(' ',1)
    Write (messages(1),'(1x,a)') 'If the above parallelization settings deteriorate the performance of the simulation,&
                                & the user might consider'
    If ((.Not. hpc_data%cpus_per_node%fread) .And. (.Not. hpc_data%number_nodes%fread)) Then
      Write (messages(2),'(1x,a)') 'specifying directives "number_nodes" and "cpus_per_node" depending on the requested&
                                  & total numbers of CPUs and the HPC architecture.'
    Else If ((.Not. hpc_data%cpus_per_node%fread) .And. hpc_data%number_nodes%fread) Then
      Write (messages(2),'(1x,a)') 'specifying directives "cpus_per_node" and adjusting "number_nodes" depending on the&
                                  & requested total numbers of CPUs and the HPC architecture.'
    Else If (hpc_data%cpus_per_node%fread .And. hpc_data%number_nodes%fread) Then
      Write (messages(2),'(1x,a)') 'adjusting directives "cpus_per_node" and "number_nodes" depending on the&
                                  & requested total numbers of CPUs and the HPC architecture.'
    End If
    Call info(messages,2)
    
    If (Trim(hpc_data%parallelism_type%type) == 'mpi-openmp') Then
      Write (messages(1),'(1x,a)') 'The user should also optimise the value of "threads_per_process".'
      Call info(messages,1)
    End If

    If (Trim(hpc_data%platform%type) == 'slurm') Then
      Write (messages(1),'(1x,2a)') 'For more information about SLURM, please check ', Trim(web_slurm)
      Call info(messages,1)
    End If

  End Subroutine summary_hpc_settings

  Subroutine get_hpc_time(input_time, timeformat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Generate time limit from the input values 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: input_time(3) 
    Character(Len=256), Intent(  Out) :: timeformat

    Integer(Kind=wi)   :: i
    Character(Len=256) :: chartime(3), aux
 
    Do i= 2, 3
      If (input_time(i)<10) Then
        Write (aux,*) input_time(i)
        chartime(i)='0'//Adjustl(Trim(aux))
      Else
        Write (aux,*) input_time(i)
        chartime(i)=Adjustl(Trim(aux))
      End If 
    End Do 

    If (input_time(1)>0) Then
      Write (aux,*) input_time(1)
      chartime(1)=Adjustl(Trim(aux))
      timeformat=Trim(chartime(1))//'-'//Trim(chartime(2))//':'//Trim(chartime(3))//':00'
    Else
      timeformat=Trim(chartime(2))//':'//Trim(chartime(3))//':00'
    End If
 
  End Subroutine get_hpc_time

End Module hpc_setup 
