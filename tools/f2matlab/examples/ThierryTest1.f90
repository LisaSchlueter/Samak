module gblcom
  ! record format of output files of routine oftrigger
   type :: trigger_record
     sequence
     integer*8 :: stamp         ! time stamp [us since epoch]
     integer   :: det           ! detector
     integer   :: fnum          ! file number
     real      :: ph            ! pulse height from optimal filter
     real      :: ofrms         ! rms from optimal filter 
     real      :: random_amp    ! amplitude of superimposed std event (at random time)
     integer   :: imax          ! array index of maximum of shaped pulse
  end type trigger_record


  ! record format of output files of routine merge_stamp_files
  type :: merge_record
     sequence
     integer*8 :: gstamp
     integer   :: det
     integer   :: file_number
     real      :: ph
     real      :: ofrms
     real      :: random_amp
     real      :: tpa
     integer   :: tpch
     integer   :: tp_delay
     integer   :: trigger_delay
     real      :: dead_time     
  end type merge_record
  
  !***************** parameters ***********************************************
  integer, parameter :: &
       KHDM=10, &              ! integers in event header        (dim. ih)
       KADM=8, &               ! dvm data in event header	 (dim. dvm)
       KTDM=4, &               ! number of times in header       (dim. tt)
       NDETM=64, &             ! number of detectors
       NDIGCHM=8, &            ! number of input channels per digitizer module   
       KPA=100, &              ! number of parameters 	         (dim. ip)
       num_flags=640, &        ! number of available flags
       MAXDIGCH=8,&            ! number of channels per digitizer module
       stream_only=1, &        ! data structure (only streaming with torso .par file))
       lux_ccs=2, &            ! data structure of CRESST DAQ at Gran Sasso
       muc_ccs=3, &            ! Munich Linux based DAQ
       KPR=360, &              ! number of standard events and filters
       maxargc=50,&            ! largest number of parameters
       argv_maxlen=60 ,&         ! max length of argv string

      ! *** parameter definitions moved from cmp to gblcom ********
       nptime=1,&        ! time since start [h]
       npcoinc=2,&       ! number of coincident pulses
       nptpa=3,&         ! amplitude of test pulses
       npfile=4,&        ! file number
       npnev=5,&         ! event number
       nptotdelay=6,&    ! delay of trigger to first trigger in event
       npheater=7,&      ! output of control DAC [V] = heating level
       nptpdelay=8,&     ! delay of test pulse trigger to fired tp (ms)
       npqdc_tim=9,&     ! time since muon [us]
       npqdc_sum=10,&    ! muon pulse height
       npqdc_mult=11,&   ! muon multiplicity
       nptpch=12,&       ! test pulse channel which fired tpa
       npcm =18, &       ! num of parameters common in event       

       npdet=npcm+1,&
       nplive=npcm+2,&
       npdelay=npcm+3,&
       npbase=npcm+4,&
       npph=npcm+7,&
       nponset=npcm+8,&
       npofamp=33,&
       npofrms=34,&
       npoframp=kpa-1 
  
  real, parameter :: nonsense=1.e36
  !****************************************************************************
  REAL :: &
       DVM(KADM),&                  ! dvm readings in event header
       TT(KTDM),&                   ! times in event header
       TIME_BASE,&                  ! time base in micro secs
       tcalfact=1.,&                ! factor for calibrating time stamping clock
       stamp_time_base,&            ! time base of time stamping clock 
       tsum,&                       ! total measuring hours in all files in list
       parfile_version              ! version number stated in .par file

  real, dimension(ndetm) :: &
       stab_timu,&         ! trigger delay time limit in stab cut 
       stab_timo,&         ! trigger delay time limit in stab cut
       stab_limu,&         ! lower limit of cntrl pulse height for stab. cut
       stab_limo,&         ! upper limit of cntrl pulse height for stab. cut
       stab_live_removed,& ! live hours removed by stability cut
       cvset,&             ! control set point
       dsum                ! total dead hours in all files in list
       

  real, dimension(kpa,ndetm) :: &
       gauge, &                      ! gauge factor
       xyl, xyh, &                   ! display section for parameter(detector)
       xyl_complete, xyh_complete, & ! complete section
       xyl_selected, xyh_selected    ! selected section


  integer, dimension(ndetm) :: &
       KAM, &                 ! channels in moving average window
       koff_limit, &          ! shift base line calculation in cmp up to this limit
       kofi, &                ! offset (baseline) calculation from in CMP
       koff, &                ! offset (baseline) calculation to in CMP
       kmaf, &                ! calculate maximum up to this channel in cmp
       triggers               ! number of triggers in digitizer channel

   integer :: &
       ih(khdm), &         ! integer data in header
       ndet, &                ! largest detector number found in cmp
       npc, &                 ! number of parameters common calculated by CMP
       npd, &                 ! number of parameters per detector calculated CMP
       krd, &                 ! actual record length of digitizer
       kra, &                 ! maximal reading of digitizer 
       nlfl, &  
       nast, & 
       nevtat, &              ! totoa number of events
       rectat, &              ! total number of record
       khd, &                 ! number of integer data in event header (includes ulong)
       khd_ulong, &           ! unsigned longs in in event header
       kad, &                 ! number of dvm data in event header
       ktd, &                 ! number of time entries in event header 
       mode, &                !
       pre_trigger, &
       len_namin, &           ! number of characters in namin
       len_namint, &          ! number of characters in namint
       len_namdir, &          ! number of chars in namdir (directory of input files)
       num_files, &           ! number of input files in file list
       tra_channels, &        ! number of active digitizer channels in tra
       tra_bits, &
       trigger_mode, &        ! 0=separate, 1=or mode, 2=group mode
       data_structure         ! lux_ccs, muc_ccs,  vms_ccs   

  logical :: &
       control_pulses_saved, & ! wave form of control pulses save in .rdt file
       sevs_modified, &        ! sevs may have to be saved on exit
       stamps_modified,&       ! time stamps may have to be saved on exit
       use_streaming_data,&    ! read streaming data instead of records from conv. DAQ
       use_ofstamps            ! use stamps produced by oftrigger (opt filter trigger)



  CHARACTER :: &
       PTX(KPA)*30="", &                ! text for paremeters
       FTX(0:num_flags)*30="", &        ! text for flags
       STX(kpr)*30="", &                ! text for standard events
       HTX(kpr)*30="", &                ! text for shaping filters
       sdet*60, &                       ! string for reading detector names
       sdetx*60, &                      !   "
       sdety*60, &                      !   "
       sdetz*60, &                      !   "
       NAMIN*88, &                      ! name of infile including directory no extendion 
       NAMH*10, & 
       namint*44, &                     ! list file, directory removed
       namdir*44, &                     ! directory of list file
       progdir*64,&                     ! directory which contains rop
       savedir*80,&                     ! directory for saving parameters, stamps and flags
       streamdir*40,&                   ! directory containing streaming data
       stream_prepend*20,&              !
       titstr(0:ndetm)*80,&             ! detector specific display titles read from par file
       stab_ignore_single(ndetm)*1      ! ignore isolated excursion flag in stability cut


  !************************ allocatable arrays *********************************
  integer*8, allocatable :: stamp(:)   ! time stamp clock values (restart at file begin)
  integer*8, allocatable :: &
       t_file_begin(:),t_file_end(:),&   ! unix time at begin and end of file
       stream_delay(:),&                 ! start delay of streaming DAQ for file
       stab_usecs_file(:,:)              ! measurement time [us] removed by stability cut  
  integer*4, allocatable :: flag(:,:)

  real*4, allocatable :: &
       dead_time_to_add(:,:),&
       stab_dead_secs_file(:,:),&        ! dead secs removed by stability cut (ndetm,num_files)
       time_to_add(:),&
       cntrl_tintv(:)


  integer*4, allocatable :: &
       recsum(:),&
       cntrl_pulses_per_testpulse(:)

  real*4, allocatable :: &             ! SR(KRD,KPR),&               ! Standard events
       sr(:,:)
  
  character*80, allocatable :: flist(:)
  

  external fip, ip, stpi, stpf ! memory managment for parameters storeparam.c
  real*4 fip
  integer*4 ip


contains
  subroutine alloc_num_files()
    !*******************************************************************
    !allocate memory for arrays with bounds depending on number of 
    !files num_files
    !********************************************************************
    implicit none
    integer alloc_status
    allocate(t_file_begin(num_files), t_file_end(num_files),&
         stream_delay(num_files),&
         time_to_add(num_files),&
         cntrl_tintv(num_files),&
         recsum(0:num_files),&
         cntrl_pulses_per_testpulse(num_files),&
         flist(num_files),&

         dead_time_to_add(ndetm, num_files),&
         stab_usecs_file(ndetm,num_files),&
         stab_dead_secs_file(ndetm,num_files),&
         stat=alloc_status)
     if(alloc_status.ne.0) then
       print *,"Error allocating memory for arrays depending on number of files"
       stop
    endif
    stream_delay(:)=0_8
    time_to_add(:)=0.
    cntrl_tintv(:)=0
    t_file_begin(:)=0_8
    t_file_end(:)=0_8
    recsum(:)=0
    cntrl_pulses_per_testpulse=0
    flist(:)=""

    dead_time_to_add(:,:)=0
    stab_usecs_file(:,:)=0
    stab_dead_secs_file(:,:)=0
  end subroutine alloc_num_files


  subroutine allocate_stamps()
    implicit none
    integer alloc_status
    allocate(stamp(rectat), stat=alloc_status)   ! memory for time stamps
    if(alloc_status.ne.0) then
       print *,"Error allocating memory for time stamps"
       stop
    endif
  end subroutine allocate_stamps

  subroutine allocate_std_events()
    implicit none
    integer alloc_status
    allocate(sr(krd,kpr), stat=alloc_status)
    if(alloc_status.ne.0) then
       print *,"Error allocating memory for standard events"
       stop
    endif
  end subroutine allocate_std_events
end module gblcom
