module utils
  character, parameter :: esc*1=char(27)
  character, parameter :: bold_green*9=esc//"[1m"//esc//"[32m"
  character, parameter :: green*5=esc//"[32m"
  character, parameter :: bold_yellow*9=esc//"[1m"//esc//"[33m"
  character, parameter :: yellow*5=esc//"[33m"
  character, parameter :: bold_magenta*9=esc//"[1m"//esc//"[35m"
  character, parameter :: magenta*5=esc//"[35m"
  character, parameter :: bold_red*9=esc//"[1m"//esc//"[31m"
  character, parameter :: red*5=esc//"[31m"   
  character, parameter :: white*5=esc//"[37m"
  character, parameter :: black*5=esc//"[30m"
  !  character, parameter :: reset*4=esc//"[0m"

  !reset bg black and fg white and reset bold
  character, parameter, private :: reset*14=esc//"[0m"//esc//"[40m"//esc//"[37m"
  character(len=*),parameter,private :: supported_color(10)=&
       (/bold_green, green, bold_yellow, yellow, bold_magenta, magenta,&
       bold_red, red, white, black/)
  external IsFlagSet
  logical IsFlagSet
contains

  subroutine set_background_and_forground_colors()
    print '(a,$)',reset
  end subroutine set_background_and_forground_colors

  subroutine reset_colors()
    print'(a,$)',esc//"[0m"
  end subroutine reset_colors

  !**********************************************************************
  !                      cprint
  !**********************************************************************
  subroutine cprint(text, which_color, new_line)
    ! print text bracketed by # in color, # are not printed
    implicit none
    character(len=*), intent(in) :: text
    character(len=*), optional, intent(in) :: which_color
    logical, optional, intent(in) :: new_line
    integer :: i, n
    logical col, nl
    character*10 color
    color=bold_green
    if(present(which_color)) then
       do i=1,size(supported_color)
          if(which_color.eq.supported_color(i)) then
             color=supported_color(i)
             exit
          endif
       end do
    endif

    nl=.true.
    if(present(new_line)) nl=new_line ! Zeilenvorschub nach Textausgabe

    col=.false.
    n=0
    do i=1,len_trim(text)
       if(text(i:i)=='#') then
          col=.not.col
          if(col) then
             print '(a,$)',trim(color)
          else
             print '(a,$)',reset
          endif
       else
          if(n>0) then
             print '(a,$)',text(i:i)
          else
             print '(1x,a,$)',text(i:i)
          endif
          n=n+1
       endif
    end do
    if(nl) print *
  end subroutine cprint

  subroutine cnprint(txt)
    character(len=*) txt
    call cprint(txt,green,.false.)
  end subroutine cnprint

  subroutine ceprint(txt)
    character(len=*), intent(in) :: txt
    call cprint(txt,bold_red,.true.)
  end subroutine ceprint


  subroutine cenprint(txt)
    character(len=*), intent(in) :: txt
    call cprint(txt,bold_red,.false.)
  end subroutine cenprint

  !********************************************************************
  !                   PrintHelp
  !********************************************************************
  subroutine PrintHelp(hlpfile)
    use gblcom
    implicit none
    character*(*) hlpfile
    CHARACTER ANS*1, txt*96, color*10
    integer ios, lun
    ans=" "
    open(newunit=lun,file=trim(progdir)//trim(hlpfile),iostat=ios)
    if(ios.ne.0) then
       print *,"Error opening file "//trim(progdir)//trim(hlpfile)
       return
    endif
    do
       read(lun,"(A)",iostat=ios) txt
       if(index(txt,"*****").gt.0) then
          color=bold_magenta
       elseif(index(txt,"-----").gt.0) then
          color=bold_yellow
       else
          color=bold_green
       endif
       if(ios.ne.0) exit
       call cprint(trim(txt),color)
    end do
    close(unit=lun)
    return
  end subroutine PrintHelp

  !***************************************************************************
  !     aski, askf, askd, aska fragen mit dem 'text' nach der Eingabe.
  !     Der vorgeschlagene Wert bei <cr> beibehalten.
  !***************************************************************************
  subroutine aski(text,nvar)
    implicit none
    CHARACTER(len=*), intent(in) :: text
    integer nvar
    integer len, ios
    character*20 an*20
    do
       write(an,'(i10)') nvar
       call cnprint(trim(text)//"  "//trim(an)//" ?  ")

       read(*,'(a)',iostat=ios) an
       len=Len_Trim(an)
       if(len > 0) then
          read(an(1:len),*,iostat=ios) nvar
          if(ios.ne.0) then
             call ceprint('#illegal input#')
             cycle
          endif
       endif
       exit
    end do
    return
  end subroutine aski

  subroutine askf(text,xvar)
    implicit none
    character(len=*), intent(in) :: text
    character*20 an
    real xvar
    integer len,ios
    do
       write(an,"(f12.4)") xvar
       call cnprint(trim(text)//" "//trim(an)//" ?  ")
       read(*,'(a)',iostat=ios) an
       len=len_trim(an)
       if(len.gt.0) then
          read(an(1:len),*,iostat=ios) xvar
          if(ios.ne.0) then
             call ceprint("#WRONG INPUT !!!!#")
             cycle
          endif
       end if
       exit
    end do
    return
  end subroutine askf

  subroutine askd(text,xdvar)
    implicit none
    character*(*) text
    real*8 xdvar
    real xvar
    xvar=sngl(xdvar)
    call askf(text,xvar)
    xdvar=xvar
  end subroutine askd

  subroutine aska(text,avar)
    implicit none
    character(len=*), intent(in) :: text
    character(len=*) avar
    integer len, ios
    character*80 an
    call cnprint(trim(text)//"  "//trim(avar)//" ?  ")
    read(*,"(a)",iostat=ios)  an
    len=len_trim(an)
    if(len>0) then
       call MakeUpperCase(an)
       avar=an
    endif
    return
  end subroutine aska


  subroutine ask_text(text,avar,upper)
    ! upper selects whether answer is converted to upper case
    implicit none
    character(len=*), intent(in) :: text
    character(len=*) avar
    logical, intent(in) :: upper
    integer len, ios
    character*80 an
    call cnprint(trim(text)//"  "//trim(avar)//" ?  ")
    read(*,'(a)',iostat=ios)  an
    len=len_trim(an)
    if(len>0) then
       if(upper) call MakeUpperCase(an)
       avar=an
    endif
    return
  end subroutine ask_text

  subroutine askc(avar)
    ! just read string and convert to upper case
    ! do not propose any present value of string
    implicit none
    character(len=*) avar
    integer ios,len
    character*80 an
    read(*,'(a)',iostat=ios) an
    len=len_trim(an)
    if(len>0) then
       call MakeUpperCase(an)
       avar=an
    endif
    return
  end subroutine askc

  !*****************************************************************************
  !                            asd
  !****************************************************************************
  logical function asd(text, detstr, det)
    !asd(text, detstr, det) 
    !asd returns true when valid detector selected
    !asd0 returns also true when ch=0 is selected and does
    !not display error message
    USE gblcom
    implicit none
    character(len=*), intent(in) ::text
    character(len=*), intent(inout) :: detstr
    integer, intent(inout) :: det
    integer ios, k, n
    character(len=len(sdet)) str, prtstr(ndetm)

    logical asd0
    logical skip_warning
    skip_warning=.false.
1   asd=.false.
    str=""
    if(skip_warning) then     ! ch0 is allowed
       print "(1x,a,$)",trim(text)//" (0=all,999=info) ["//&
            trim(detstr)//"]: "
    else
       print "(1x,a,$)",trim(text)//" (search string, number or 999 for list) ["//trim(detstr)//"]: "
    endif


    read(*,"(a)",iostat=ios) str
    if(len_trim(str).eq.0) str=detstr

    read(str,"(i4)",iostat=ios) det 
    if(ios.eq.0) then         ! input is number
       if(det.eq.999) then
          do k=1,ndetm
             if(triggers(k).eq.0) cycle
             if(len_trim(titstr(k)).eq.0) cycle
             print "(i4,a)",k,": "//trim(titstr(k))
          end do
          goto 1
       elseif(det.gt.ndetm) then 
          print *,char(7)," Detector number to large."
          print *," Largest possible detector number is ",ndetm
          return
       endif
       if(det.lt.1) then
          if(.not. skip_warning) then
             print *,char(7)," Detector number too small."
             print *," Smallest possible number  is",ndetm
             return
          else
             det=0
             asd=.true.
             return
          endif
       endif
       asd=.true.
       detstr=titstr(det)
    else                      ! input is string
       n=0
       do k=1,ndetm
          if(triggers(k).eq.0) cycle
          if(len_trim(titstr(k)).eq.0) cycle
          if(len_trim(str).eq.0) cycle
          if(iindex(titstr(k),trim(str)).eq.0) cycle
          det=k
          n=n+1
          prtstr(n)=titstr(k)(1:len(prtstr(n)))
       end do
       if(n.eq.0) then 
          print *," String "//trim(str)//&
               " does not match any detector name"
          goto 1
       endif
       if(n.eq.1) then
          detstr=prtstr(1)
          asd=.true.
          return
       endif
       if(n.gt.1) then
          do k=1,n
             print "(a,i3,a)","Channel",k,", Detector name="//trim(prtstr(k))
          end do
          print *,"More than one detector name includes search string."//&
               " Enter unique search string or number."
          goto 1
       endif
    endif
    return

    entry asd0(text, detstr, det)
    skip_warning=.true.
    goto 1 
  end function asd



  !*********************************************************************
  !                           askh
  !**********************************************************************
  subroutine askh(nh, h_write)
    !ask which transfer function to use
    !input:     
    !nh       suggested value for transfer function
    !h_write  when set will askh will ask for new text
    !output:    
    !nh  accepted value. nh=0 indicates failure (nothing was accepted)

    USE gblcom
    implicit none
    integer nh
    logical  h_write
    integer ios
    character*1 ans
    do
       call aski("Number of transfer function (999=Info):",nh)
       if (nast.gt.0) then
          nh=0
          return
       endif

       if(nh.eq.999) then
          call ptl("")
          cycle
       endif

       if(nh.lt.1 .or. nh.gt.kpr) then
          print *,"Number outside accepted range of 1 to ",kpr
          cycle
       endif
       if(.not.h_write) then
          if(len_trim(htx(nh)).eq.0) then
             print *,"Selected transfer function is not defined"
             cycle
          endif
          print "(a,/)"," Selected transfer function: "//trim(htx(nh))
          return
       endif
       exit
    end do

    !now we have a valid number for a transfer function to be overwritten

    if(len_trim(htx(nh)).gt.0) then
       ans="N"
       print *," Transfer function already defined: "//trim(htx(nh))
       call aska("Overwrite this transfer function: ",ans)
       if(ans.eq."N") then
          nh=0
          return
       endif

       ans="N"
       call aska("Delete text :",ans)
       if(ans.ne."Y") return 
    endif


    do
       print *,"Enter text for transfer function"
       print *,"<--------------------------->"
       read (5,"(a)",iostat=ios) htx(nh)
       if(ios.ne.0 .or. len_trim(htx(nh)) .eq. 0) then
          print *,"Error reading text"
          cycle
       endif
       exit
    end do


    return
  end subroutine askh


  !****************************************************************
  !                         asp
  !****************************************************************
  SUBROUTINE ASP(TEXT,NP)
    ! ASP(TEXT,NP) FRAGT MIT "TEXT" NACH PARAMETERNUMMER NP
    USE gblcom
    implicit none
    integer np
    CHARACTER*(*) TEXT  
    integer len_text, len_npstr, ios
    character npstr*30, ans*1
    len_text=index(text,":")
    if(len_text.gt.0) then
       len_text=len_text-1
    else
       len_text=len_trim(text)
    endif
101 if(nast.gt.0) return
    print "(a,i3,a,$)",text(1:len_text)//" (text or 999=info) [",np,"]: "
    read (*,"(a)",iostat=ios) npstr
    if(nast.gt.0) return

    len_npstr=len_trim(npstr)
    if(len_npstr.eq.0) then
       if(np.gt.kpa .or. np.lt.1) then
          print "(a,i3,a)", "Parameter number ",np," out of range"
          print "(a,i1,a,i3,a)","Numbers from ",1," to ",kpa," are allowed"
          goto 101
       elseif(len_trim(ptx(np)).gt.0) then
          return
       else
          print "(a,i3,a)","Parameter number ",np," is empty"
          ans="Y"
          call aska("Would you like to use this empty parameter:",ans)
          if(nast.gt.0) return
          if(ans.ne."Y") goto 101
          return
       endif
    endif
    read(npstr,*,iostat=ios) np
    if(ios.ne.0) then  ! np contains text
       call ppl(npstr)
       goto 101
    endif
    IF(NP.EQ.999) THEN
       CALL PPL("")
       GOTO 101
    endif

    if(np.gt.kpa .or. np.lt.1) then
       PRINT *," Only 1 to ",KPA," permitted !"
       GOTO 101                                                    
    ENDIF
    if(len_trim(ptx(np)).eq.0) then
       print "(a,i3,a)","Parameter number ",np," is empty"
       ans="Y"
       call aska("Would you like to use this empty parameter:",ans)
       if(ans.ne."Y") goto 101
    endif
    RETURN
  END SUBROUTINE ASP

  !*************************************************************************
  !                       ASF
  !*************************************************************************
  SUBROUTINE ASF(TEXT,IC,NFL)
    ! ASF(TEXT,IC,NFL) fragt nach flag nummmer.
    ! IC=0: LESEFLAG 
    ! IC=1: SETZFLAG
    ! IC=2: SETZFL. UEBERSCREIBEN ohne Rückfrage
    ! IC=-1:LESEFLAG -1 MOEGLICH
    ! FLAGTEXTE WERDEN AUSGEGEBEN BZW. EINGELESEN.
    USE gblcom
    implicit none
    integer ic, nfl, ifl, ios, n
    character*(*)text
    integer len_text
    character nflstr*30, ans*1

    len_text=index(text,":")
    if(len_text.gt.0) then
       len_text=len_text-1
    else
       len_text=len_trim(text)
    endif

401 call cprint("Enter #999# to list all flags or #a search string#")
    print "(a)","to list flags with text containing this string"
    if(nfl.lt.10) then
       print "(a,i1,a,$)", text(1:len_text)//" [",nfl,"]: "
    elseif(nfl.lt.100) then
       print "(a,i2,a,$)", text(1:len_text)//" [",nfl,"]: "
    elseif(nfl.lt.1000) then
       print "(a,i3,a,$)", text(1:len_text)//" [",nfl,"]: "
    else
       print "(a,i4,a,$)", text(1:len_text)//" [",nfl,"]: "
    endif

    read(*,"(a)",iostat=ios) nflstr
    if(nast.gt.0) return

    if(ios.ne.0) call ceprint("#Error reading nflstr#")
    if(len_trim(nflstr).eq.0) then
       if(nfl.eq.999) then
          call pfl("")
          print *
          goto 401
       endif
    else
       read(nflstr,*,iostat=ios) nfl
       if(ios.ne.0) then !nflstr may contain flag text to print selected list
          call pfl(nflstr)
          print *
          goto 401
       endif
    endif

    IF(NFL.EQ.999) THEN! print complete flag list
       CALL PFL("")
       print *
       GOTO 401
    ENDIF

    IF(NFL.GT.num_flags.OR.(NFL.LT.0.AND.IC.NE.-1)) GOTO 401
    IF(NFL.LT.0) GOTO 999
    if(ic.eq.0 .and. len_trim(ftx(nfl)).eq.0.) then
       call ceprint("#This flag is empty !!!!#")
       goto 401
    endif
    if(ic.eq.1 .or. ic.eq.2) then
       if(nfl.eq.0) then
          call ceprint("#Overwriting of flag 0 not possible#")
          goto 401
       endif
    endif

    CALL GFL(NFL,IFL)
    IF(IC.LE.0.OR.IFL.NE.0.OR. len_trim(FTX(NFL)).gt.0) THEN
       PRINT 291, NFL,IFL,FTX(NFL)
291    FORMAT(/,1X,"Flag",I4," contains ",1X,I9," events.",/,1X,"Text: ",A30,/)
       IF(IC.EQ.2) THEN
          PRINT *,"Will be deleted !"
          GOTO 505
       ENDIF
    ENDIF
    IF(NFL.EQ.0.OR.IC.LE.0) return

    IF(IFL.gt.0) then
       ANS=" "
       CALL ASKA("Flags (#CR#=)Keep / (#D#)elete :",ANS)
       IF(ANS.EQ."D") THEN
          DO N=1,NEVTAT
             CALL SFL(N,NFL,0)
          end do
       ENDIF
    endif
    IF(len_trim(FTX(NFL)).gt.0) then
       ANS=" "
       CALL ASKA("Text  (#CR#=)Keep / (#D#)elete :",ANS)
       IF(ANS.ne."D") return
    endif
505 PRINT 510, NFL
510 FORMAT(/,"<----- Text for Flag ",I3," ---->")
    read (*, "(A30)") FTX(NFL)
999 RETURN                           
  END SUBROUTINE ASF



  !**********************************************************************
  !                          PFL
  !**********************************************************************
  SUBROUTINE PFL(name)
    !************************************************************
    ! pfl prints the number of events with flags set
    ! When name is a string with nonzero length only flags
    ! with a title are listed which contains this string
    !*************************************************************
    use gblcom
    implicit none
    character*(*) name
    integer :: n, ifl, nn, nn_was, len_name, nfound

    nn_was=0
    nfound=0
    len_name=len_trim(name)
    if(len_name.eq.0) then
       PRINT 105
       PRINT 36, 0,nevtat,FTX(0)
    endif
    DO N=1,num_flags 
       if(len_trim(ftx(n)).eq.0) cycle
       if(nast.eq.1) exit
       if(len_name.gt.0) then
          if(iindex(ftx(n),name(1:len_name)).gt. 0) then
             call gfl(n,ifl)
             nn=n/10
             if(nn.ne.nn_was) then
                if(nfound.gt.0)&
                     print *,"******************************************"
                nn_was=nn
             endif
             PRINT 36, N,IFL,FTX(N)
             nfound=nfound+1
          endif
       else
          CALL GFL(N,IFL)
          if(ifl.ne.0) then
             nn=n/10
             if(nn.ne.nn_was) then
                print *,"******************************************"
                nn_was=nn
             endif
             PRINT 36, N,IFL,FTX(N)
          endif
       end if
    end do
105 FORMAT(1X,"Flag  Content    <---------- Text ------------>")
36  FORMAT(1X,I4,"  ",I8,"   ",A30)
    if(len_name.gt.0 .and. nfound.eq.0) then
       print *," There is no flag text which contains this string: "//&
            name(1:len_name)
    endif
    RETURN
  END SUBROUTINE PFL

  !**************************************************************
  !           GET_FLAG_LIST asks for a list of flags
  !**************************************************************
  subroutine get_flag_list(text, nfl, n, nmax)
    !---------------------------------------------------------
    !  input:
    !    text   text to be printed
    !    nfl    the flag list (will also be output when question
    !           is prompted with only return)
    !    n      number of flags that was obtained (not altered when
    !           question prompted with return)
    !    nmax   length of nfl in calling program
    !
    !  output:
    !    nfl, n
    !----------------------------------------------------------
    use gblcom
    implicit none
    integer n, nmax
    integer nfl(nmax)
    character*(*) text
    integer ios, i, len
    character*80 nflstr, ans
    integer nfl_temp(nmax)

    write(nflstr,*,iostat=ios)(nfl(i),", ",i=1,n)
    if(ios.ne.0) then
       call ceprint("#Error writing flag list to string#")
    endif
    call compact_string(nflstr)
    len=max(1,len_trim(nflstr)-1)
401 print "(a,$)",text(1:len_trim(text))//" ["//nflstr(1:len)//"]: "
    read(*,"(a)", iostat=ios) ans
    len=len_trim(ans)
    if(len.gt.0) then
       read(ans,*,iostat=ios)(nfl_temp(i),i=1,n)
       if(.not. GetListedNumbers(ans, nfl_temp, n, nmax, 0, 999,.false.)) then
          call pfl(ans)
          goto 401
       end if

       do i=1,n
          if(nfl_temp(i).eq. 999) then
             call pfl("")
             goto 401
          endif
          if(nfl_temp(i) .gt. num_flags) then 
             print *," Flag number larger than allowed maximum "
             goto 401
          endif
          if(nfl_temp(i) .lt. 0) then
             print *," Flag number < 0 not allowed"
             goto 401
          endif
       end do
       do i=1,n
          nfl(i)=nfl_temp(i)
       end do
    else
       do i=1,n
          if(nfl(i).eq. 999) then
             call pfl("")
             goto 401
          endif
          if(nfl(i) .gt. num_flags) then 
             print *," Flag number larger than allowed maximum "
             goto 401
          endif
          if(nfl(i) .lt. 0) then
             print *," Flag number < 0 not allowed"
             goto 401
          endif

       end do

    endif
    return
  end subroutine get_flag_list

  !**********************************************************************
  !                  IsInList
  !**********************************************************************
  logical function IsInList(n,list,length)
    implicit none
    integer n,length
    integer list(length)
    integer i
    do i=1,length
       if(n.eq.list(i)) then
          IsInList=.true.
          return
       endif
    end do
    IsInList=.false.
    return
  end function IsInList

  !**********************************************************************
  !                      IsAnyFlagSet
  !**********************************************************************
  logical function IsAnyFlagSet(nev,nfl,n)
    implicit none
    integer nev, n
    integer nfl(n)
    integer i
    external IsFlagSet
    logical IsFlagSet
    IsAnyFlagSet=.true.
    do i=1,n
       if(IsFlagSet(nev,nfl(i))) return
    end do
    IsAnyFlagSet=.false.
    return
  end function IsAnyFlagSet




  !**********************************************************************
  !                     fipg
  !**********************************************************************
  real function fipg(rec,npar)  
    USE gblcom, only : fip, ip, npdet, gauge
    implicit none
    integer, intent(in) :: rec, npar
    fipg=fip(rec,npar)*gauge(npar,ip(rec,npdet))
    return
  end function fipg

  !**********************************************************************
  !                       mktime
  !**********************************************************************
  subroutine mktime(systime, time_string)
    ! convert system time in printable askii string containing UTC time
    integer*8 systime
    character*(*) time_string
    integer t(9)
    integer*8 system_time_8
    integer system_time(2)
    equivalence(system_time(1),system_time_8)
    character*3 day(0:6) /'Sun', 'Mon','Tue','Wed','Thu','Fri','Sat'/
    character*3 month(0:11) /'Jan','Feb','Mar','Apr','May','Jun','Jul',&
         'Aug','Sep','Oct','Nov','Dec'/

    system_time_8=systime !this defines system_time
    if(system_time(1).gt.1000000) then
       print *,' low byte of time stamp > 1000000'
    endif
    call gmtime(system_time(2), t)
    write(time_string,1000) day(t(7)), month(t(5)), t(4), t(3), t(2), t(1),&
         system_time(1)/1000, 1900+t(6)
1000 format(a3,1x,a3,1x,i2.2,1x,i2.2,':',i2.2,':',i2.2,'.',i3.3,' UTC ',i4.4)
    return
  end subroutine mktime
  !**********************************************************************
  !                            get_unix_time
  !**********************************************************************
  integer*8 function get_unix_time(tim)
    !     convert time in us into Unix format which contains seconds in high 
    !     word and us in low word
    implicit none
    integer*8 tim             ! time in [us]
    integer*8 t
    integer*4 tw(2)
    integer*8, parameter :: mio=1000000

    equivalence(t,tw)

    tw(2)=int(tim/mio)
    tw(1)=int(tim - tw(2)*mio)
    get_unix_time=t
    return
  end function get_unix_time




  !**********************************************************************
  !                            GetCommands
  !**********************************************************************
  integer function GetCommands(argv, argcmax, script, nast)
    !*************************************************************************
    !tokenize command line
    !input: 
    !     argcmax      the largest number of tokens allowed by calling program
    !output:
    !      argv        Array of Null terminated tokens. argv[0] is converted
    !                  to upper case
    !The function returns the number of argv tokens obtained
    !*************************************************************************

    implicit none
    integer, intent(in) :: argcmax,  nast
    character(len=*), intent(in) :: argv(argcmax)
    character(len=*), intent(inout) :: script    !script to run at start of rop
    integer len, argc

    character(len=240)  in    
    external getlin
    integer getlin

    !get_lin is called with len set to size of string and returns actual lenth
    GetCommands=0
    argc=0
    len=getlin(in,trim(script)//char(0),nast)
    script="";
    if(len.eq.0) return
    call strtok(in,argv(1))
    if(argv(1)(1:1).eq.' ') return
    argc=1
    do
       if(argc.eq.argcmax) then
          call cprint("#Too many elements in argument list. List truncated#")
          exit
       endif
       call strtok(char(0),argv(argc+1))
       if(argv(argc+1)(1:1).eq.' ') exit
       argc=argc+1
    end do
    GetCommands=argc
    return
  end function GetCommands


  !**********************************************************************
  !                        get_option_value
  !**********************************************************************

  integer function get_option_value(option, value, argc, argv)
    !
    !     GetOptionValue 
    !     Purpose: Extract value for option from any of the argv fields
    !     
    !     input :
    !        option   The option string to search for. If option is "-d" then a
    !                 string containing part of a detector name may follow. I this
    !                 case the returned value is the detector name that contains
    !                 the strig following -d
    !        arg!     number of argv elements
    !        argv     string containing option and value
    !     output: 
    !        value    The value which follows option, It can be either in the same 
    !                 argv element or in the element following the option string
    !     
    !     A non zero return value indicates successful extraction of value
    !******************************************************************************
    use gblcom
    implicit none
    integer rget_option_value
    character(len=*) option
    real value
    integer argc
    character(len=*) argv(argc)
    logical :: rem
    integer i, k, lopt, largv, ios
    character*20 bargv

    rem=.false.
1   lopt=Len_Trim(option)
    get_option_value=0
    i=0
    do while (i.lt.argc)
       i=i+1
       !options start with -
       if(argv(i)(1:lopt+1).ne.'-'//trim(option)) cycle    
       largv=Len_Trim(argv(i))
       if(largv .gt. lopt+1) then !value in same token
          read(argv(i)(lopt+2:largv),*,iostat=ios) value
          if(ios.eq.0) then
             if(rem) call remove_the_argument(i,1)
             get_option_value=1 ! success
             return
          endif
          if(option(1:1).eq.'d') then ! might be detector name
             bargv=argv(i)(2+lopt:largv)
             !              print *,' searching title strings for '//trim(bargv)
             do k=1,ndetm
                if(iindex(titstr(k),trim(bargv)).gt.0) then
                   value=k
                   get_option_value=1
                   if(rem) call remove_the_argument(i,1)
                   return
                endif
             end do
          endif
       else                   !value in next token
          if(argc.gt.i) then  ! next token available
             i=i+1
             read(argv(i),*,iostat=ios) value
             if(ios.eq.0) then
                get_option_value=1 ! success
                if(rem) call remove_the_argument(i-1,2)
                return
             endif
             if(option(1:1).eq.'d') then ! might be detectorr name
                largv=len_trim(argv(i))
                do k=1,ndetm
                   if(iindex(titstr(k),trim(argv(i))).gt.0) then
                      value=k
                      get_option_value=1
                      if(rem) call remove_the_argument(i-1,2)
                      return
                   endif
                end do
             endif
          endif
       endif
    end do
    return
    entry  rget_option_value(option, value, argc, argv)
    rem=.true.
    goto 1
  contains 
    subroutine remove_the_argument(nfrom,nremove)
      integer nfrom,nremove
      argc=argc-nremove
      do k=nfrom,argc
         argv(k)=argv(k+nremove)
      end do
      return
    end subroutine remove_the_argument
  end function get_option_value


  !**********************************************************************
  !                   get_option_string
  !**********************************************************************
  integer function get_option_string(option, value, argc, argv)
    !Purpose: Extract string value for option from any of the argv
    ! Input:
    !  option   The option string to search for, initial - stripped
    !    argc   number of argv elements
    !    argv   string containing option and value
    !  Output:
    !    value string which follows option. The value string can be either in 
    !          the same or in the following argv element.
    !A non zero return value indicates successful extraction of value
    use gblcom
    implicit none
    integer, intent(in) :: argc
    character(len=*), intent(in) :: option, argv(argc)
    character(len=*), intent(out) :: value
    integer i, k, lopt, largv

    lopt=Len_Trim(option)
    get_option_string=0
    i=0
    do while (i.lt.argc)
       i=i+1
       if(argv(i)(1:1+lopt) .ne. '-'//trim(option)) cycle                 ! options has to start with -
       largv=Len_Trim(argv(i))
       if(largv .gt. lopt+1) then                    !value in same token
          value=argv(i)(lopt+2:largv)
          get_option_string=1 ! success        
          do k=i+1,argc                              !first token starting with - terminates the string
             if(argv(k)(1:1).eq.'-') return
             value=trim(value)//' '//argv(k)
          end do
          return
       else                                          !value in next token
          do k=i+1,argc
             if(k.eq.i+1) then
                !              if(argv(k)(1:1).eq.'-') then
                !                 get_option_string=0 ! error
                !                 print *,'Error: Option instead of text is following '//&
                !                      '-'//trim(option)
                !                 return
                !              endif
                value=argv(k)
                get_option_string=1 ! success
             else
                if(argv(k)(1:1).eq.'-')  return
                value=trim(value)//' '//argv(k)
             endif
          end do
       endif
    end do
    return
  end function get_option_string

  !**********************************************************************
  !                         options_valid
  !**********************************************************************
  logical function options_valid(nopts, option,argc,argv)
    !check if all command line options are valid
    !input:
    !  nopts         number of supported options
    !  option        supported options with leading - stripped
    !  argc. argv    command line arguments
    implicit none
    integer, intent(in) :: nopts, argc
    character(len=*), intent(in) :: option(nopts), argv(argc)
    integer i, n, lena, leno
    logical valid
    options_valid=.true.
    do n=1,argc
       if(argv(n)(1:1).ne.'-') cycle
       lena=max(last_alpha(argv(n)),index(argv(n),"("), index(argv(n),")"))
       valid=.false.
       do i=1, nopts
          leno=len_trim(option(i))
          if(lena.eq.leno+1) then
             if(option(i)(1:leno).eq.argv(n)(2:lena)) then
                valid=.true.
                exit
             endif
          endif
       end do
       if(.not.valid) then
          options_valid=.false.
          call ceprint("#Option "//trim(argv(n))//" not supported#")

       end if
    end do
    return
  end function options_valid

  integer function last_alpha(txt)
    ! returns position of last alpha (A-Z, a-z) in string or 0 in case of
    ! failure
    ! skips possibly leading -- or -
    implicit none
    character(len=*), intent(in) :: txt
    integer n, nfrom, len

    last_alpha=0
    len=len_trim(txt)
    if(len.lt.2) return
    nfrom=1
    if(txt(1:2).eq."--") then
       nfrom=3
    else if(txt(1:1).eq."-") then
       nfrom=2
    endif

    do n=nfrom,len
       if(isalpha(txt(n:n))) then
          last_alpha=n
       else 
          exit
       endif
    end do
  end function last_alpha



  !***********************************************************
  !                      option_present
  !***********************************************************
  integer function option_present(option, argc, argv)
    ! !!!!! enter option with leading - stripped
    implicit none
    integer, intent(in) :: argc
    character(len=*) option, argv(argc)
    integer i, lopt, largv
    option_present=1
    do i=1, argc
       largv=last_alpha(argv(i))
       lopt=len_trim(option)
       if(largv.eq.lopt+1) then
          if(argv(i)(1:largv).eq."-"//option(1:lopt)) return
       endif
    end do
    option_present=0
    return
  end function option_present


  !**********************************************************************
  !                         isalpha
  !**********************************************************************
  logical function isalpha(c)
    ! retruns true if c is a-z or A-Z, false otherwise
    implicit none
    character*1, intent(in) :: c
    integer n
    n=ichar(c)
    isalpha=.false.
    if((n.ge.65 .and. n.le.90) .or. ((n.ge.97 .and. n.le.122))) isalpha=.true.
  end function isalpha

  logical function isnumber(c)
    ! retruns true if c is 0-9, false otherwise
    implicit none
    character*1, intent(in) :: c
    integer n
    n=ichar(c)
    isnumber=.false.
    if((n.ge.48) .and. (n.le.57)) isnumber=.true.
  end function isnumber



  !**********************************************************************
  !                           get_option_list
  !**********************************************************************
  integer function get_option_list(option, value, argc, argv, nel, nelmax)
    !Purpose: Extract list of values for option from any of the argv fields     
    !input :
    !  option   The option string to search for. 
    !  argc     number of argv elements
    !  argv     string containing option and value
    !  nelmax
    !output:
    !  nel      Number of values returned
    !  value    The value which follows option, It can be either in
    !           the same argv element or in the following element 
    !Returns number of values obtained
    implicit none
    integer, intent(in) :: argc, nelmax
    character(len=*), intent(in) :: option, argv(argc)
    real, intent(out) :: value(nelmax)
    integer, intent(out) :: nel 

    integer :: i, k, lopt, largv, ios, dashid
    real val1, val2
    character str*80
    get_option_list=0
    nel=0
    lopt=Len_Trim(option)
    i=0
    do
       i=i+1
       if(i.gt.argc) exit
       if(argv(i)(1:lopt+1).ne. '-'//option(1:lopt)) cycle 
       str=argv(i)(2:)                          
       largv=len_trim(str)
       if(largv .gt. lopt) then
          str=str(lopt+1:)       
          largv=len_trim(str)  
       else
          i=i+1
          str=argv(i)
       endif


       do                             ! get list of values
          dashid=index(str,'-')
          if(dashid.gt.1) then
             read(str(1:dashid-1),*,iostat=ios) val1
             if(ios.ne.0) exit
             read(str(dashid+1:),*,iostat=ios) val2
             if(ios.ne.0) exit
             do k=int(val1), int(val2)
                nel=nel+1
                value(nel)=k
             end do
          else                         ! values are in following tokens
             read(str,*,iostat=ios) val1
             if(ios.ne.0) exit
             nel=nel+1
             value(nel)=val1
          endif
          i=i+1
          if(i.gt.argc) exit
          str=argv(i)                          
       end do
       exit
    end do
    get_option_list=nel
    return
  end function get_option_list


  !******************************************************************************
  !                            GetListedNumbers
  !******************************************************************************
  logical function GetListedNumbers(list_string, list, n, nmax, &
       listval_min, listval_max, sorted)
    !-------------------------------------------------------------------------
    ! input:
    !     list_string    string from which listed numbers shall be extracted
    !     nmax           length of array list in calling program
    !     listval_min    allowed minumum of values in list
    !     listval_max    allowed maximum of vulues in list
    !     sorted         specifies if values in list are sorted (default=.true. 
    ! output:
    !     list           integer array with resulting list
    !     n              number of elements in list
    !
    ! Rhe function returns .true. on success and .false. on error
    !
    ! Accepted formats of list string: 
    !   1-11,20    without spaces
    !   -11,20     get all values from listval_min to -11 and 20
    !   11-,20     get all values from 11 to listval_max and 20
    !
    !list is sorted and and list elements occurring more than once are removed
    !-------------------------------------------------------------------------
    implicit none
    character(len=*), intent(in) :: list_string
    integer, intent(in) :: nmax, listval_min, listval_max
    integer, intent(out) :: list(nmax), n
    logical, optional, intent(in) :: sorted
    integer i, k, l, dash_pos, kfrom ,kto, list_min, kmin,&
         listval_from, listval_to

    integer sorted_list(nmax)     ! for sorting
    logical sorting
    logical mask(nmax)            ! for sorting
    character(len=len_trim(list_string)) str
    sorting=.true.
    if(present(sorted)) sorting=sorted

    GetListedNumbers=.false.
    if(len_trim(list_string).eq.0) return


    if(list_string.eq.'*') then
       GetListedNumbers=.true.
       k=listval_min
       do i=1,nmax
          list(i)=k
          if(k.eq.listval_max) then
             n=k
             print *,' all elements selected. List length is ',n
             return
          endif
          k=k+1
       end do
    endif

    i=1
    kfrom=1
    do 
       str=list_string(kfrom:)
       kto=index(str,',')
       if(kto.eq.0) kto=len_trim(str) ! no further comma present
       if(kto.eq.0) exit     

       dash_pos=index(str(1:kto),'-')
       if(dash_pos.gt.0) then
          if(dash_pos.eq.1) then      
             read(str(2:),*,end=999,err=999) listval_to
             listval_from=listval_min
          elseif(dash_pos.eq.kto) then      !dash at end of last token
             read(str(1:kto-1),*,end=999,err=999) listval_from
             listval_to=listval_max 
          elseif(dash_pos.eq.kto-1) then    
             if(str(kto:kto).eq.',') then
                read(list_string(1:kto-2),*,end=999,err=999) listval_from
                listval_to=listval_max                
             else
                read(list_string(1:dash_pos-1),*,err=999,end=999) listval_from
                read(str(dash_pos+1:),*,err=999,end=999) listval_to
             endif
          else
             read(list_string(1:dash_pos-1),*,err=999,end=999) listval_from
             read(str(dash_pos+1:),*,err=999,end=999) listval_to
          end if

          do 
             list(i)=listval_from
             i=i+1
             if(i.gt.nmax) then
                print *,' Too many elements in list'
                return
             endif
             if(listval_from.eq.listval_to) exit
             listval_from=listval_from+1
          end do

       else
          read(str,*,end=999, err=999) list(i)
          i=i+1
          if(i.gt.nmax) then
             print *,' Too many elements in list'
             return
          endif
       endif
       kfrom=kfrom+kto
    end do

    n=i-1
    GetListedNumbers=.true.
    do i=1,n
       if(list(i).gt. listval_max .or. list(i) .lt. listval_min) then
          print *,' Some value in list of Numbers is out of range'
          GetListedNumbers=.false.
          exit
       end if
    end do

    !check for double occurence of list elements and remove

    i=1
    do while(i.lt.n)
       k=i+1
       do while(k.le.n)
          if(list(i).eq.list(k)) then
             n=n-1
             do l=i,n
                list(l)=list(l+1)
             end do
          else
             k=k+1
          endif
       end do
       i=i+1
    end do

    if(.not.sorting) return

    !sort in ascending order


    mask(1:n)=.true.


    do i=1,n
       list_min=2**30
       kmin=0
       do k=1,n
          if(mask(k).and.(list(k).le.list_min)) then
             kmin=k
             list_min=list(k)
          endif
       end do
       if(kmin.gt.0) then
          mask(kmin)=.false.
          sorted_list(i)=list(kmin)
       endif
    end do
    do i=1,n
       list(i)=sorted_list(i)
    end do
    return
999 print *,' Error reading list of numbers'
    return
  end function GetListedNumbers


  !**********************************************************************
  !*                  remove_multiple_spaces (replace them by single)
  !**********************************************************************
  subroutine remove_multiple_spaces(str)
    character(len=*), intent(inout) :: str
    integer slen
    integer i,k
    slen=len(str)
    k=1
    do i=1,slen-1
       if(str(i:i+1).eq."  ") cycle
       str(k:k)=str(i:i)
       k=k+1
    end do
    do i=k,slen-1
       str(i:i)=" "
    end do
    return
  end subroutine remove_multiple_spaces

  !********************************************************************
  !                   compact_string   (remove spaces)
  !********************************************************************
  subroutine compact_string(str)
    character(len=*) str
    integer i,k,slen
    slen=len(str)
    k=1
    do i=1,slen
       if(str(i:i).eq.' ') cycle
       str(k:k)=str(i:i)
       k=k+1
    end do
    do i=k,slen
       str(i:i)=char(32)
    end do
    return
  end subroutine compact_string

  !**********************************************************************
  !                    MakeUpperCase
  !**********************************************************************      
  subroutine MakeUpperCase(str)
    implicit none
    character*(*) str
    integer i,k,n
    k=len_trim(str)
    do i=1,k                 ! convert command to upper case
       n=ichar(str(i:i))
       if(n.ge.97.and.n.le.122) str(i:i)=char(n-32)
    end do
    return
  end subroutine MakeUpperCase
  !**********************************************************************
  !                     UpperCase
  !**********************************************************************
  function UpperCase(str)
    !trim string and make it upper case
    implicit none
    character(Len=*), intent(in) :: str
    character(len=len_trim(str)) :: UpperCase
    integer i,n,k

    k=len_trim(str)
    UpperCase=str(1:k)
    do i=1,k
       n=ichar(str(i:i))
       if(n.ge.97.and.n.le.122) UpperCase(i:i)=char(n-32)
    end do
    return
  end function UpperCase

  !**********************************************************************
  !                     count_tokens
  !**********************************************************************
  subroutine count_tokens(instr,count)
    !count comma or space separated tokens
    implicit none
    character*(*), intent(in) :: instr
    integer, intent(out) :: count
    character(len=len_trim(instr)) str
    integer len, i, n

    count=0
    str=adjustl(instr)
    len=len_trim(str)
    if(len.eq.0) return

    do i=1,len-1
       if(str(i:i).eq.",") str(i:i)=" "
    end do
    n=0
    do i=1,len-1
       if((str(i:i).eq." ").and.(str(i+1:i+1) .ne. " ")) n=n+1
    end do
    count=n+1
    return
  end subroutine count_tokens

  !*******************************************************************
  !  iindex is case insensitive version of intrinsic function index
  !*******************************************************************
  integer function iindex(str, substr)
    ! find position of trimmed string substr in str
    implicit none
    character(len=*), intent(in) :: str, substr
    iindex=index(UpperCase(str), UpperCase(substr))
    return
  end function iindex


  !**********************************************************************
  !             us_to_cont_hours
  !**********************************************************************
  real function us_to_cont_hours(time_stamp, filenum)
    use gblcom, only : t_file_begin, t_file_end
    implicit none
    integer*8, intent(in) :: time_stamp
    integer, intent(in) :: filenum
    integer*8 tsum, t
    integer n
    tsum=0
    if(filenum>1)  then
       n=filenum-1
       tsum=sum(t_file_end(1:n) -t_file_begin(1:n))
    endif
    t = (time_stamp - t_file_begin(filenum)  + tsum)/10000
    us_to_cont_hours=float(t) / 3600.e2
  end function us_to_cont_hours



  !**********************************************************************
  !             convert_clock_to_us
  !**********************************************************************
  subroutine convert_clock_to_us(stamp, file_number)
    ! return time stamp entered in clock units as time stamp in
    ! us since epock
    use gblcom, only : t_file_begin, stamp_time_base
    implicit none
    integer*8, intent(inout) :: stamp
    integer, intent(in) :: file_number
    if(stamp_time_base .eq. 0.1) then
       stamp=stamp/10_8
    elseif(stamp_time_base.eq.10.) then
       stamp=stamp*10_8
    else
       call cenprint("#Value of stamp_time_base is not supported: #")
       print *,stamp_time_base
       stop
    endif
    stamp=stamp+t_file_begin(file_number)
  end subroutine convert_clock_to_us

  !**********************************************************************
  !                    cloc_units_to_us
  !**********************************************************************
  integer*8 function clock_units_to_us(time_stamp,filenum)
    use gblcom, only : stamp_time_base, t_file_begin
    implicit none
    integer*8, intent(in) :: time_stamp
    integer, intent(in) :: filenum

    if(stamp_time_base .lt. 1) then
       clock_units_to_us = time_stamp/10 + t_file_begin(filenum)
    else 
       clock_units_to_us=time_stamp*10 + t_file_begin(filenum)
    endif
  end function clock_units_to_us

  !**********************************************************************
  !                     us_to_clock_units
  !**********************************************************************
  integer*8 function us_to_clock_units(time_stamp, filenum)
    use gblcom, only : stamp_time_base, t_file_begin
    implicit none
    integer*8, intent(in) :: time_stamp
    integer, intent(in) :: filenum

    if(stamp_time_base .lt. 1) then
       us_to_clock_units=(time_stamp-t_file_begin(filenum))*10
    else
       us_to_clock_units=(time_stamp-t_file_begin(filenum))/10
    endif
  end function us_to_clock_units


  !**********************************************************************
  !                    rtd
  !**********************************************************************
  subroutine rtd(argc, argv)
    !    read transfer functions from disk
    use gblcom, only : krd, kpr, htx, savedir
    implicit none
    integer, intent(in) :: argc
    character(len=*), intent(in) :: argv(argc)
    integer i, ifrom, ito, n, ios, nh, k, lun, lun1, lenh, saved_nh
    logical rkeep(kpr)
    character*80 fname
    character*42 temp_file(2) 


    if(ShowHelp(argc,argv)) then
       call cprint("\n&
            &******************************* RTD *******************************************\n&
            &RTD serves to read transfer functions. Without options, specific transfer\n&
            &functions in the results directory (files ftf_nnn) are read after those from\n&
            &files generic.ftf_nnn in the working directory. Functions imported from generic\n&
            &files are overwritten by specific ones of the same number. Some options are\n& 
            &available to change this default behaviour\n&
            & Sntax: RTD [options]\n&
            &#  Options:#\n&
            &#  -g    #Read filters only from files generic.ftf_nnn and skip specific\n&
            &#  -s    #Read only specific filters from result directory\n&
            &#  -k    #Keep (do not overwrite) existing transfer functions of the present\n&
            &#        #rop session\n&
            &*******************************************************************************")
       return
    endif


    temp_file(1)='temp_directory_listing_1.txt'
    temp_file(2)='temp_directory_listing_2.txt'

    ifrom=1
    ito=2
    if(option_present("g",argc,argv).eq.1) then
       ito=1
       if(option_present("s",argc,argv).eq.1) then
          call ceprint("#Option -s ignored due to presence of option -g#")
       endif
    else if(option_present("s",argc,argv).eq.1) then
       ifrom=2
    end if

    rkeep(:)=.false.
    if(option_present("k",argc,argv).eq.1) then
       do k=1,kpr
          rkeep(k)=(len_trim(htx(k)).gt.0)
       end do
    else
       htx(1:kpr)=""
    endif

    if(ifrom.eq.1) call system("ls -1 generic.ftf_* > "//trim(temp_file(1))//" 2>/dev/null")  
    if(ito.eq.2) call system('ls -1 '//trim(savedir)//'/ftf_* > '//trim(temp_file(2))//" 2>/dev/null")       
    do i=ifrom,ito
       open(newunit=lun,file=temp_file(i),status='old',iostat=ios)
       if(ios.ne.0) then
          call ceprint("#Error opening file "//trim(temp_file(i))//"#")
          cycle
       endif
       do k=1,kpr
          read(lun,'(a)',iostat=ios) fname
          if(index(fname,"~").gt.0) cycle
          if(ios.ne.0) exit
          nh=get_filter_number(fname)
          if(nh.le.0) cycle
          if(rkeep(nh))then
             print *,'keeping '//htx(nh)
             cycle
          endif

          open(newunit=lun1,file=fname,form='unformatted',access='stream',&
               iostat=ios, status='old')
          if(ios.ne.0) then
             call ceprint("#File "//trim(fname)//" does not exist#")
             return
          endif


          read(lun1,iostat=ios) lenh, saved_nh
          if(lenh .ne. 2*krd) then
             print *,' Filter in file '//trim(fname)//' has incorrect length: ',lenh
          else
             read(lun1,iostat=ios) n,htx(nh)(1:n-1)
             if(ios.ne.0) then
                print *,'Error reading filter from file '//trim(fname)
                htx(nh)=""
             else
                print '(a,i3,a)'," Filter nr.",nh," read from "//trim(fname)
             endif
          endif
          close(unit=lun1)

       end do
       close(unit=lun)
       call unlink(temp_file(i),ios)
       if(ios.lt.0) then
          call perror("failure deleting file "//trim(temp_file(i)))
       endif
    end do
    return
  contains
    !**********************************************************************
    !                   
    !**********************************************************************
    integer function get_filter_number(infile)
      ! returns filter number found or -1 on failure
      implicit none
      character (len=*), intent(in) :: infile
      integer nh, n, ios
      get_filter_number=-1

      nh=-1
      n=index(infile,'_',back=.true.)
      if(n.gt.0) then
         read(fname(n+1:),*,iostat=ios) nh
         if(ios.ne.0) then
            print *,'Error extracting appended number string of file '//trim(fname)
            return
         endif
      else
         call ceprint("#No appended number string in file "//trim(fname)//"#")
         return
      endif
      get_filter_number=nh
    end function get_filter_number
  end subroutine rtd

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !                         write_filter
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  logical function write_filter(nh, H, nsev, rav, sav, filter_type)
    use gblcom, only : savedir, htx, krd
    implicit none
    integer, intent(in) :: &
         
         nh, &         ! filter number
         nsev          ! number of std event from which it was built
    real, intent(in) :: &
         H(0:2*krd-1),&  ! transfer function of filter
         rav(0:krd-1),&  ! normalized template from std event
         sav(0:krd/2)    ! power spectrum of noise
    character*1, intent(in) :: filter_type

    character*80 fname
    integer ios, n, lun
    write_filter=.false.

    write(fname,"(a,i3.3)") trim(savedir)//"/ftf_",nh
    open(newunit=lun,file=fname,form="unformatted",access="stream",iostat=ios)
    if(ios.ne.0) then
       call ceprint("#File "//trim(fname)//" can not be opened#")
       return
    endif
    n=len_trim(htx(nh))
    write(lun,iostat=ios) &
         2*krd,&                ! length of filter
         nh,&                   ! rop specific number of filter
         n+1,&                  ! length of text describing filter including "\0" 
         htx(nh)(1:n),char(0),& ! text describing filter
         h(0:2*krd-1),&         ! the filter
         nsev,&                 ! number of std event it was built from
         rav(0:krd-1),&         ! normalized std event    (not zero extended)  
         sav(0:krd/2),&         ! power spectrum of noise (not zero extended)
         filter_type
    if(ios.ne.0) then
       call ceprint("#Error occurred while writing to "//trim(fname)//"#")
       goto 999
    endif
    print *,"Filter saved in file: "//trim(fname)
    write_filter=.true.
999 close(unit=lun)
  end function write_filter


  !******************************************************************
  !              wsd: write standard events to disk
  !******************************************************************
  subroutine wsd(argc, argv)
    !default is saving in specific file only, when some std events
    !have change
    use gblcom
    implicit none
    integer, intent(in), optional :: argc
    character(len=*), optional :: argv(*)
    integer k, n, ndef, num_def, ios, lun
    logical file_exists, with_arguments
    character fname*80
    character*20 version_string /"STD Version 1.0"/
    with_arguments=(present(argc) .and. present(argv))
    if(with_arguments) then
       if(ShowHelp(argc,argv)) then
          call PrintHelp("wsd.hlp")
          return
       endif
    endif
    if(len_trim(ptx(5)).eq.0) then
       call ceprint("#Run CMP or RPD before you think about saving#")
       return
    endif

    ndef=0
    do n=1,kpr
       if(len_trim(stx(n)).gt.0) ndef=ndef+1
    end do
    if(ndef.eq.0) then
       call ceprint("#No std events defined. Nothing to save#")
       return
    endif

    fname=trim(savedir)//"/std_events"
    if(with_arguments) then
       if(option_present("f",argc,argv) .eq. 1) then
          sevs_modified=.true.
       else if(option_present("g",argc,argv) .eq.1) then
          sevs_modified=.true.
          fname="generic.std"
       endif
    else if(.not. sevs_modified) then
       print *,'Std. events did not change and will not be saved'
       return
    endif

    inquire(file=fname, exist=file_exists)
    if(file_exists) then
       if(.not.MakeBackupCopy(trim(fname))) then
          call ceprint("#Failure making backup copies of "//trim(fname)//"#")
       endif
    endif
    open(newunit=lun,file=fname,form="UNFORMATTED", status="NEW",iostat=ios)
    if(ios.ne.0) then
       call ceprint("#Error opening file "//trim(fname)//"#")
       call ceprint("#Standard events can not be saved#")
       return
    endif

    write(lun,iostat=ios) version_string
    num_def=0
    do k=1,kpr
       if(len_trim(stx(k)).gt.0) then
          write(lun,iostat=ios) k
          write(lun,iostat=ios) stx(k)
          write(lun,iostat=ios) sr(1:krd,k)
          num_def=num_def+1
       endif
    end do
    close(unit=lun)
    if(ios.eq.0) then
       print "(a,i3)","Number of standard events written to "//&
            trim(fname)//": ", num_def
       sevs_modified=.false.
    else
       call ceprint("#Error writing std. events to file "//trim(fname)//"#")
    endif
    return
  end subroutine WSD

  !******************************************************************
  !             RSD: Read Standard event form Disk
  !******************************************************************
  subroutine RSD (argc, argv)
    ! read first std events from generic files and then add/overwrite
    ! with those found in specific file
    use gblcom
    implicit none
    integer, intent(in) :: argc
    character(len=*), intent(in) :: argv(argc)
    real version_number
    real dummy(krd)
    integer i, k, kfrom, kto, nrd, num_def, lun, ios
    logical skeep(kpr)
    character(len=len(stx)) txt
    character*80 fname(2)*80, version_string_rd*20

    if(ShowHelp(argc,argv)) then
       call PrintHelp("rsd.hlp")
       return
    endif

    fname(1)="generic.std"
    fname(2)=trim(savedir)//"/std_events"
    kfrom=1
    kto=2

    if(option_present("g",argc,argv).eq.1) then
       kto=1
       if(option_present("s",argc,argv).eq.1) then
          call ceprint("#Option -s ignored due to presence of the -g option#")
       endif
    else if(option_present("s",argc,argv).eq.1) then
       kfrom=2
    endif

    skeep(:)=.false.
    if(option_present("k",argc,argv).eq.1) then
       do i=1,kpr
          if(len_trim(stx(i)).gt.0) then
             skeep(i)=.true.
          else
             stx(i)=""
             sr(1:krd,i)=0.
          endif
       end do
    else
       stx(1:kpr)=""
       sr(1:krd,1:kpr)=0.
    endif

    num_def=0
    do k=kfrom,kto
       open(newunit=lun,file=fname(k),STATUS='OLD',FORM='UNFORMATTED', iostat=ios)
       if(ios.ne.0) cycle

       read(lun,iostat=ios) version_string_rd
       if(ios.ne.0) then
          print *," Error reading version string from file "//trim(fname(k))
          goto 999
       endif

       version_number=0.
       if(index(version_string_rd,"STD Version").ne.0) then
          read(version_string_rd(12:),*,iostat=ios) version_number
       endif

       if(version_number.eq. 0.) then
          rewind(lun)
          nrd=60 ! this was kpr before introduction of version string
          if(kpr.lt.nrd) then
             print *,"kpr (max nr. of std events) too small to read std events"
             print *,"Present value of kpr: ",kpr
             print *,"Std events in file (nrd): ", nrd
             print *,"Compile with kpr>=nrd to make file readable"
             goto 999
          endif

          stx(nrd:kpr)=""
          read(lun,iostat=ios)stx(1:nrd)
          if(ios.ne.0) then
             print *,"Error reading text fields for standard events"
             print *,"File: "//trim(fname(k))
             goto 999
          endif
          num_def=0
          do i=1,kpr
             if(len_trim(stx(i)).gt.0) then
                read(lun,iostat=ios) sr(1:krd,i)
                if(ios.ne.0) then
                   call ceprint("#Error reading standard events from disk#")
                   call ceprint("#File: "//trim(fname(k))//"#")
                   goto 999
                endif
                num_def=num_def+1
             endif
          end do
       else
          do i=1,kpr
             read(lun,iostat=ios) nrd
             if(ios.ne.0) exit
             if(nrd.lt.1 .or. nrd.gt.kpr) then
                print *," Std event number out of range: ",nrd
                ios=-1
                exit
             endif
             if(skeep(nrd)) then
                ! make dummy read
                print *,"kept existing std event: "//trim(stx(nrd))
                read(lun,iostat=ios) txt 
                if(ios.ne.0) exit
                read(lun,iostat=ios) dummy
                if(ios.ne.0) exit

             else
                read(lun,iostat=ios) stx(nrd)
                if(ios.ne.0) exit
                read(lun,iostat=ios) sr(1:krd,nrd)
                if(ios.ne.0) exit
                num_def=num_def+1
             end if
          end do
       endif
    end do
    print '(a,i4)'," Nr. of standard events read: ",num_def
    sevs_modified=.false.
999 close(unit=lun)
  end subroutine RSD

  !**************************************************************************
  !     MakeBakupCopy Make Backup copy when file exists
  !**************************************************************************
  logical function MakeBackupCopy(outfile)
    implicit none
    character(len=*), intent(in) :: outfile
    integer, parameter :: num_backup=2
    integer len_outfile, len_fn, k, ios
    logical file_exists
    character(len=len_trim(outfile)+num_backup) :: fn
    len_outfile=len_trim(outfile)
    MakeBackupCopy=.false.

    fn=outfile(1:len_outfile)
    len_fn=len_outfile
    do k=1,num_backup
       fn=fn(1:len_fn)//'~'
       len_fn=len_fn+1
    end do
    !    print *

    do k=1,num_backup
       inquire(file=fn(1:len_fn-1), exist=file_exists)
       if(file_exists) then
          call rename(fn(1:len_fn-1),fn(1:len_fn),ios)
          if(ios.ne.0) then
             call ceprint("#Creation of "//fn(1:len_fn)//" failed#")
          else
             MakeBackupCopy=.true.
          endif
       endif
       len_fn=len_fn-1
    end do
  end function MakeBackupCopy



  !**********************************************************************
  !                          asks
  !**********************************************************************

  subroutine asks(nsev,sev_write, optional_text)
    !ask which standard event to select
    !input: 
    !       nsev   suggested number of standard event
    !       sev_write  when set select std event for writing
    !output: accepted value of standard event
    !       nsev=0 indicates that nothing was accepted

    USE gblcom
    integer, intent(inout) :: nsev
    logical, intent(in) :: sev_write
    character(len=*), optional, intent(in) :: optional_text

    integer ios, len_text
    logical available
    character str*30, text*80, ans*1

    if(present(optional_text)) then
       text=optional_text
    else
       text='Nr. of standard-event (text or 999=info):'
    endif
    len_text=len_trim(text)

    if(.not.sev_write) then      
       available=.false.
       do i=1,kpr
          if(len_trim(stx(i)).ne.0) available=.true.
       end do
       if(.not.available) then
          call ceprint("#There is no standard event available#")
          call ceprint("#Use the command SSE (Sum Standard Event) to define one#")
          call ceprint("#or the command RTD to read available ones from disk#")
          nsev=0
          return
       endif
    endif
    !  call psl                          ! list available std events

20  print '(a,i4,a,$)',text(1:len_text)//' ',nsev,' ? '
    read(*,'(a)',iostat=ios) str
    if(nast.ne.0) return
    if(ios.ne.0) then
       print *,' Error reading from stdin'
       goto 20
    endif
    len_str=len_trim(str)
    if(len_str.gt.0) then
       read(str,*,iostat=ios) nsev
       if(ios.ne.0) then      ! input is string
          call psl(str)
          goto 20
       endif
    endif

    if(nsev.eq.999) then
       call psl('')
       goto 20
    endif

    IF(NSEV.GT.KPR) then
       call cenprint("#Largest possible number:#")
       print "(i5)",kpr
       goto 20
    endif

    if(nsev.lt.1) then
       call cenprint("#Smallest possible number: 1#")
       GOTO 20
    endif

    if(sev_write) then
       ans='Y'
       if(len_trim(stx(nsev)).gt.0) then
          call ceprint("#Standard Event already defined:# "//trim(stx(nsev)))
          ans='N'
          call aska(' Overwrite Standard Event: ',ans)
          if(ans.eq.'Y') then
             ans='N'
             call aska(' Delete text :',ans)
          else  
             nsev=0
             return
          endif
       endif
       if(ans.eq.'Y') then
          ios=1
          do while (ios.ne.0)
             print '(a)','<------- Text--------------->'
             read (*,'(a)',iostat=ios) stx(nsev)
             if(ios.ne.0) print *,char(7),'Error reading text'
          end do
       endif
    else
       if(len_trim(stx(nsev)).eq.0) then                   
          print *,char(7),' Standard event not defined'
          goto 20        
       endif
       print *,'selected standard event: '//stx(nsev)
       print *
    endif
  end subroutine asks

  !*******************************************************
  !                askst
  !*******************************************************
  subroutine askst(optional_text,nsev,sev_write)
    ! like asks but with optional text for asking
    implicit none
    character(len=*), intent(in) :: optional_text
    integer, intent(inout) :: nsev
    logical, optional, intent(in) :: sev_write
    logical swrite
    if(present(sev_write)) swrite=sev_write
    call asks(nsev,swrite,optional_text)
  end subroutine askst

  !********************************************
  !             print standard event list
  !********************************************
  subroutine psl(name)
    USE gblcom
    implicit none
    character*(*) name
    integer len_name, n, k
    len_name=len_trim(name)
    print *,' **************************************************'
    n=0
    do k=1,kpr
       if(len_name.gt.0) then
          if(iindex(stx(k),name(1:len_name)).ne.0) then
             n=n+1
             print '(2x,i7,5x,a)',k,stx(k)
          endif
       else
          if(len_trim(stx(k)).gt.0) then
             print '(2x,i7,5x,a)',k,stx(k)
          ENDIF
       endif
    end do
    print *,' **************************************************'
    if(len_name.gt.0) then
       if(n.eq.0) then
          print *,char(7)//'There is no std event which contains this text: '//&
               name(1:len_name)
       endif
    else
       print '(A,I4)',' Largest number of standard events: ',kpr
    endif
    return
  end subroutine psl


  !**************************************************************
  !                      ppl
  !**************************************************************
  SUBROUTINE PPL(text)
    !PPL DRUCKT DIE TEXTE ZU DEN PARAMETRN AUS
    USE gblcom
    implicit none
    character*(*) text
    integer len_text, k, n

    len_text=len_trim(text)
    if(len_text .eq. 0) PRINT 643
    n=0
    DO K=1,KPA
       if(len_text.ne.0) then
          if(iindex(ptx(k),text(1:len_text)).gt.0)then
             PRINT 642, K,PTX(K)
             n=n+1
          endif
          cycle
       endif
       if(len_trim(ptx(k)).gt.0) then
          PRINT 642, K,PTX(K)
       endif
    end do
    if(len_text.gt.0 .and. n.eq.0) then
       print *,' There is no parameter text which contains this string: '//&
            text(1:len_text)
    endif
643 FORMAT(/,1X,'Param <------------Text------------>   ')      
642 FORMAT(1X,I3,3X,A30,5X,G12.5)
    RETURN
  END SUBROUTINE PPL


  !********************************************
  !     print transfer function  list 
  !********************************************
  subroutine ptl(text)
    USE gblcom, only : kpr, htx
    character (len=*) text
    integer len_text, k
    len_text=len_trim(text)

    print *,'****** List of transfer functions ********************'
    do k=1,kpr
       if(len_trim(htx(k)).gt.0) then
          if(len_text.gt.0) then
             if(iindex(htx(k),trim(text)).gt.0) print '(2x,i7,5x,a)',k,htx(k)
          else
             print '(2x,i7,5x,a)',k,htx(k)
          endif
       endif
    end do
    print *,'**************************************************'
    print '(a,i5)',' Largest possible number of filters:',kpr
    return
  end subroutine ptl

  !**********************************************************************
  !                ShowHelp
  !**********************************************************************
  logical function ShowHelp(argc, argv)
    implicit none
    integer, intent(in) :: argc
    character*(*), intent(in) ::  argv(argc)
    ShowHelp=.false.
    if(option_present('H',argc,argv).eq.1) then
       ShowHelp=.true.
    else if(option_present('h',argc,argv).eq.1) then
       ShowHelp=.true.
    else if(option_present('-help',argc,argv).eq.1) then
       ShowHelp=.true.
    endif
  end function ShowHelp



  !**********************************************************************
  !                        sorti8
  !**********************************************************************
  subroutine sorti8(n,arr,brr)
    !sort arr and sort brr accordingly
    implicit none
    integer*8, intent(in) :: n
    integer*8, intent(inout) ::  arr(n), brr(n)

    INTEGER*8 M,NSTACK
    PARAMETER (M=7,NSTACK=50)
    INTEGER*8 i,ir,j,jstack,k,l,istack(NSTACK)
    integer*8 a,b,temp
    character*1 cont
    jstack=0
    l=1
    ir=n
1   if(ir-l.lt.M) then
       do j=l+1,ir
          a=arr(j)
          b=brr(j)
          do  i=j-1,1,-1
             if(arr(i).le.a)goto 2
             arr(i+1)=arr(i)
             brr(i+1)=brr(i)
          end do
          i=0
2         arr(i+1)=a
          brr(i+1)=b
       end do
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       temp=arr(k)
       arr(k)=arr(l+1)
       arr(l+1)=temp
       temp=brr(k)
       brr(k)=brr(l+1)
       brr(l+1)=temp
       if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
       endif
       if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
       endif
       if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
       endif
       i=l+1
       j=ir
       a=arr(l)
       b=brr(l)
3      continue
       i=i+1
       if(arr(i).lt.a)goto 3
4      continue
       j=j-1
       if(arr(j).gt.a)goto 4
       if(j.lt.i)goto 5
       temp=arr(i)
       arr(i)=arr(j)
       arr(j)=temp
       temp=brr(i)
       brr(i)=brr(j)
       brr(j)=temp
       goto 3
5      arr(l)=arr(j)
       arr(j)=a
       brr(l)=brr(j)
       brr(j)=b
       jstack=jstack+2
       if(jstack.gt.NSTACK) then
          cont='Y'
          print *,'NSTACK too small in sort2'
          call aska('Continue program (Y/N): ',cont)
          if(cont.ne.'Y') stop
       endif
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1
  end subroutine sorti8

  !**********************************************************************
  !        global     sorti81
  !**********************************************************************
  logical function sorti81(n,arr)
    INTEGER*8, intent(in) :: n
    integer*8, intent(inout) :: arr(n)
    integer*8, parameter :: M=7,NSTACK=50
    INTEGER*8 i,ir,j,jstack,k,l,istack(NSTACK)
    Integer*8 a,temp
    sorti81=.true.
    jstack=0
    l=1
    ir=n
1   if(ir-l.lt.M) then
       do j=l+1,ir
          a=arr(j)
          do i=j-1,1,-1
             if(arr(i).le.a)goto 2
             arr(i+1)=arr(i)
          end do
          i=0
2         arr(i+1)=a
       end do
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       temp=arr(k)
       arr(k)=arr(l+1)
       arr(l+1)=temp
       if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
       endif
       i=l+1
       j=ir
       a=arr(l)
3      continue
       i=i+1
       if(arr(i).lt.a) goto 3
4      continue
       j=j-1
       if(arr(j).gt.a)goto 4
       if(j.lt.i)goto 5
       temp=arr(i)
       arr(i)=arr(j)
       arr(j)=temp
       goto 3
5      arr(l)=arr(j)
       arr(j)=a
       jstack=jstack+2
       if(jstack.gt.NSTACK) then
          call ceprint("#NSTACK too small in sorti81#")
          sorti81=.false.
          return
       endif
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1
  end function sorti81

  !***********************************************************************
  !                rpf  (read parameter file)
  !**************************************************************************
  logical function rpf(fname, recs, counts, dead_time, &
       measuring_time, vset, control_pulse_interval, &
       control_pulses_per_testpulse, Id)

    use gblcom
    implicit none

    integer, intent(in) :: id
    integer, intent(out) ::  recs,  counts(ndetm), control_pulses_per_testpulse
    real, intent(out) :: dead_time(ndetm), vset(ndetm), measuring_time, &
         control_pulse_interval
    character(len=*), intent(in) :: fname


    integer*8 :: tempi8, file_size
    real z
    integer :: lun
    integer :: ios, filenum,  n, i, k, klim

    character dummy*132, numstr*20, ans*1
    logical verbose, file_exists, last_file

    verbose=.false.
    last_file=(id.eq.0) 
    rpf=.false.

    open(newunit=lun,file=trim(fname)//".par", status='old', iostat=ios)
    if(ios.ne.0) then
       call ceprint("#Parameter file not found: "//trim(fname)//".par#")
       stop
    endif

    write(dummy, "(a,i4,a)",iostat=ios) "#********************* File ",id+num_files,&
         ": "//trim(fname)//" ****************#"
    call remove_multiple_spaces(dummy)
    call cprint(trim(dummy))

    read(lun,'(a)',err=190)dummy
    if(index(dummy,"CRESST Version").eq.0) then
       call ceprint("#CRESST Version missing in first line of "//trim(fname)//".par#")
       stop
    endif
    read(dummy(16:),*) parfile_version

    stamp_time_base=0.1  
    tra_channels=ndetm        
    if(parfile_version.ge.9.0) then
       if(use_ofstamps) then
          data_structure=stream_only
       else
          data_structure=muc_ccs
          stamp_time_base=10.    ! time base of time stamping clock
          tra_channels=8        
       endif
    else
       data_structure=lux_ccs
    endif
    tra_bits=16

    !*********************************************************************
    !when only streaming data is available read only 6 lines followed by
    !detector titles (when string DISPLAY titles is present)

    ! timeofday at start [s]: 
    ! timeofday at start [us]:
    ! timeofday at stop [s]:  
    ! timeofday at stop [us]:
    ! Record length: 
    ! Time base [us]:

    ! optional string containing Display titles
    ! followed by names for each channel e.g. ch1, A-PH
    ! 
    !**********************************************************************

    if(parfile_version.gt.1.0) then
       filenum=id+num_files
       read(lun,'(a)',err=190) dummy ! timeofday at start [s]
       n=index(dummy,':')
       read(dummy(n+1:),*,err=190,end=190) tempi8

       read(lun,'(a)',err=190) dummy ! timeofday at start [us]
       n=index(dummy,':')
       read(dummy(n+1:),*,err=190,end=190) t_file_begin(filenum)
       t_file_begin(filenum)=t_file_begin(filenum)+tempi8*1000000

       read(lun,'(a)',err=190) dummy ! timeofday at stop  [s]
       n=index(dummy,':')
       read(dummy(n+1:),*,err=190,end=190) tempi8

       read(lun,'(a)',err=190) dummy ! timeofday at stop  [us]
       n=index(dummy,':')
       read(dummy(n+1:),*,err=190,end=190) t_file_end(filenum)
       t_file_end(filenum)=t_file_end(filenum)+tempi8*1000000
    endif

    if(data_structure.ne.stream_only) then
       !**************** print start writing time ***********************
       read(lun,'(a)',err=190,end=190)dummy 
       print '(1x,a,$)',trim(dummy(index(dummy,':')+1:))

       !**************** print stop writing time ************************
       read(lun,'(a)',err=190,end=190)dummy ! stop writing date and time
       print '(a)','  ===>  '//trim(dummy(index(dummy,':')+1:))


       !**************** measuring hours ********************************
       if(.not.fread_to("Measuring time",":",numstr,measuring_time)) then
          call ceprint("Error: Did not find string 'Measuring time    [h]'#")
          goto 190
       endif
       print *, " Measuring time [h]: "//trim(numstr)
       if(last_file) print *,""


       !*************** integers in event header *****************************
       if(.not.iread_to("Integers in header",":",numstr,khd)) then
          call ceprint("Error: Did not find string 'Integers in header'#")
          goto 190
       endif
       if(last_file) then
          print *," Integers in header:       "//trim(numstr)
       endif


       !*************** unsigned longs in event header *********************
       if(.not.iread_to("Unsigned longs in header",":",numstr,khd_ulong)) then
          call ceprint("Error: Did not find string 'Unsigned longs in header'#")
          goto 190
       endif
       if(last_file) then
          print *," Unsigned longs in header: "//trim(numstr)
       endif
       khd=khd+khd_ulong      ! treat ulongs as ints

       !**************** reals in event header *******************************         
       if(.not.iread_to("Reals in header",":",numstr,ktd)) then
          call ceprint("#Error: Did not find string 'Reals in header'#")
          goto 190
       endif
       if(last_file) then
          print *," Reals in header:          "//trim(numstr)! including dvm readings
       endif

       !*************** number of dvm channels *****************************
       if(.not.iread_to("DVM channels",":",numstr,kad)) then
          call ceprint("#Error: did not find string 'DVM channels'#")
          goto 190
       endif
       if(last_file) then
          print *," DVM channels:             "//trim(numstr)
       endif
       ktd=ktd-kad
    end if

    !************** record length ****************************************
    if(.not.iread_to("Record lengt",":",numstr,krd)) then
       call ceprint("#Error: did not find string 'Record lengt'#")
    endif
    if(last_file) then
       print *," Record length:            "//trim(numstr)
    endif

    if(data_structure.ne.stream_only) then
       !**************** records written ***********************************     
       if(.not.iread_to("Records written",":",numstr,recs)) then
          call ceprint("#Error: did not find string 'Records written'#")
       endif

       !***************** pretrigger **************************************
       if(.not. iread_to("Pre trigger", ":", numstr, pre_trigger)) then
          call ceprint("#Error did not find string 'Pre trigger'#")
          goto 190
       endif
       if(last_file) then
          print *," Pre trigger:              "//trim(numstr)
       endif
    endif

    !***************** time base [us]***********************************
    if(.not. iread_to("Time base",":", numstr, k)) then
       call ceprint("#Error: Did not find string 'Time base'#")
       goto 190
    endif
    time_base=k
    if(last_file) then
       print *," Time base [us]:           "//trim(numstr)
    endif



    !****************  channels groupted or paired  *******************
    if(data_structure.eq.muc_ccs) then
       if(.not.iread_to("Or trigger phonon/light",":",numstr,trigger_mode)) then
          call ceprint("#Error: Did not find string 'Or trigger phonon/light'#")
          goto 190
       endif
    elseif(data_structure.eq.lux_ccs) then
       if(.not.iread_to("Trigger mode", ":", numstr, trigger_mode)) then
          call ceprint("#Error: Did not find string 'Trigger mode'#")
          goto 190
       endif
    endif


    if(data_structure.ne.stream_only) then
       !***************** control pulses in .rdt file *********************
       read(lun,'(a)',err=190,end=190)dummy ! empty line
       read(lun,'(a)',err=190,end=190) dummy 
       control_pulses_saved=(index(dummy,'control pulses saved in separate file').eq.0)
       if(last_file) then
          if(control_pulses_saved) then
             print *," Control pulses in .rdt file"
          else
             print *," Pulse height of control pulses in separate file"
          endif
       endif
    endif


    if(.not.use_ofstamps) then   
       if(.not.fread_to("Dead time [h]","",numstr,z)) then
          call ceprint("#Error: Did not find string 'Dead time [h]'#")
          goto 190
       endif
       if(verbose) then
          print *,"***************************** dead time [%]: ****************************"
       endif
       do i=1,tra_channels/8
          read(lun,'(a)',ERR=190,END=190) dummy 
          n=index(dummy,':')
          klim=min(tra_channels-(i-1)*8,8)
          read(dummy(n+1:),*,err=190,end=190)(dead_time((i-1)*8+k),k=1,klim)   
          if(verbose) then
             print '(i2,a1,8f9.2)',(i-1)*8+1,':',&
                  (100.*dead_time((i-1)*8+k)/measuring_time,k=1,klim)
          endif
       end do
       if(verbose) then
          print *
          print *,"*****************************live time [h] ***********************************"
          do i=1,tra_channels/8
             klim=min(tra_channels-(i-1)*8,8)
             print '(i2,a1,8f9.4)',(i-1)*8+1,':',&
                  (measuring_time-dead_time((i-1)*8+k),k=1,klim)
          end do
          print *
       endif

       if(.not.iread_to("Channel counts","",numstr,k)) then
          call ceprint("Error: did not find string 'Channel counts'#")
          goto 190
       endif
       if(verbose) then
          print *,"**************************** channel counts **********************************"
       endif
       do i=1,tra_channels/8
          read(lun,'(a)',ERR=190,END=190) dummy 
          n=index(dummy,':')
          klim=min(tra_channels-(i-1)*8,8)
          read(dummy(n+1:),*,err=190,end=190)(counts((i-1)*8+k),k=1,klim)   
          if(verbose) print '(i2,a1,8i9)',(i-1)*8+1,':',(counts((i-1)*8+k),k=1,klim)
       end do
       if(verbose) print *
    endif   !.not.use_ofstamps



    if(data_structure.eq.lux_ccs .or. data_structure.eq.muc_ccs) then
       !************************ test pulses ***********************************
       if(.not.fread_to("Interval for sending test pulses [s]",":",numstr,z)) then
          call ceprint("#Error: Did not find string 'Interval for sending test pulses [s]'#")
          goto 190
       endif
       if(.not.iread_to("Control pulses per test pulse",":",numstr,i)) then
          call ceprint("#Error: Did not find string 'Control pulses per test pulse'#")
          goto 190
       endif
       control_pulses_per_testpulse=i
       control_pulse_interval=z/(i+1)

       !************ control set point [V] ***********************************
       if(.not.iread_to("Control: Set point [V]:","",numstr,i)) then
          call ceprint("#Error: did not find string 'Control: Set point [V]:'#")
          rewind(lun)
       else
          do i=1,tra_channels/8
             read(lun,'(a)',err=190,end=190) dummy
             n=index(dummy,':')+1
             klim=min(tra_channels-(i-1)*8,8)
             read(dummy(n:),*,iostat=ios) (vset((i-1)*8+k),k=1,klim) 
          end do
       endif
    endif



    !************* detector names (DAQ display titels) ********************
    if(last_file) then
       if(.not.iread_to("Display Titles" ,"",numstr,i))  then
          if(data_structure.ne.stream_only) then
             call ceprint("#Error: Did not find string 'Display Titles'#")
             goto 190
          else
             rpf=.true.
             call ceprint("#Did not find string 'Display Titles'#")
             goto 999
          endif
       endif
       titstr(0)='all detectors'
       do i=1,tra_channels
          write(titstr(i),'("ch ",I3)') i
          call remove_multiple_spaces(titstr(i))
       end do

       i=1
       ios=0
       do while (i.le.tra_channels)         
          read(lun,'(a)',iostat=ios) titstr(i)
          if(ios.ne.0) exit
          if(len_trim(titstr(i)).eq.0) exit
          if(index(titstr(i),'ch').eq.0) exit

          n=index(titstr(i),',')
          if(n.gt.0) then
             titstr(i)=adjustl(titstr(i)(n+1:))
             n=index(titstr(i),',')
             if(n.gt.1) then
                titstr(i)=titstr(i)(1:n-1)
             endif
          endif
          i=i+1
       end do
    end if

    rpf=.true.
    if(.not.use_ofstamps) then
       inquire(file=trim(fname)//".rdt",exist=file_exists, size=file_size)
       if(file_exists) then
          n=int(file_size/((KHD + KAD + KTD)*4 + KRD*2))
          if(n.ne.recs) then
             call ceprint("#!*!*!*WARNING*!*!*!* WARNING *!*!*!*!* WARNING !*!*!*!*!")
             print *, "Problem with file "//trim(fname)//".rdt"
             print *, "Number of records stated in  "//trim(fname)//".par", recs
             print *, "Number of records found in   "//trim(fname)//".rdt", n
             call ceprint("#")
             ans='Y'
             call aska(' Would you like to use value found in data file ',ans)
             if(ans.eq.'Y') then
                recs=n
             else
                stop
             endif
             call ceprint("#!*!*!*WARNING*!*!*!* WARNING *!*!*!*!* WARNING !*!*!*!*!#")
          endif
       else
          call ceprint("#File: "//trim(fname)//".rdt does not exist#")
          stop
       endif
    endif
    goto 999
190 call ceprint("#Error reading file: "//trim(fname)//".par#")
999 CLOSE(UNIT=lun)
    return

  contains

    logical function iread_to(search_string, cut_string, aval, ival)
      ! Read lines in file until line contains search_string.
      ! When cut_string is not empty aval contains sub string
      ! following cut_string and ival is integer value contained in
      ! sub string
      implicit none
      character(len=*), intent(in) :: search_string, cut_string
      integer, intent(out) :: ival
      character(len=*),intent(out) :: aval
      integer ios, i, k
      character*80 str
      iread_to=.false.
      do
         read(lun,'(a)',iostat=ios) str
         if(ios.ne.0) then
            call ceprint("iread_to: #Error: End of file encountered#")
            return
         endif
         if(index(str, trim(search_string)).gt.0) exit
      end do
      k=len_trim(cut_string)
      if(k.gt.0) then
         i=index(str,cut_string(1:k))
         read(str(i+k:),*,iostat=ios) ival
         if(ios.ne.0) then
            call ceprint("iread_to: #Error extracting int from "//trim(str)//"#")
            return
         endif

         write(aval,*,iostat=ios) ival
         if(ios.ne.0) then
            call ceprint("-read_to: #Error in internal write, writing int#")
            return
         endif
         aval=adjustl(aval)
      endif
      iread_to=.true.
    end function iread_to

    logical function fread_to(search_string, cut_string, aval, rval)
      ! Read lines in file until line contains search_string.
      ! When cut_string is not empty aval contains sub string
      ! following cut_string and ival is float in sub string
      implicit none
      character(len=*), intent(in) :: search_string, cut_string
      real, intent(out) :: rval
      character(len=*),intent(out) :: aval
      integer ios, i, k
      character*80 str
      fread_to=.false.
      do
         read(lun,'(a)',iostat=ios) str
         if(ios.ne.0) then
            call ceprint("fread_to: #Error: End of file reached#")
            return
         endif
         if(index(str, trim(search_string)).gt.0) exit
      end do
      k=len_trim(cut_string)
      if(k.gt.0) then
         i=index(str,cut_string(1:k))
         read(str(i+k:),*,iostat=ios) rval
         if(ios.ne.0) then
            call ceprint("fread_to: #Error extracting float from "//trim(str)//"#")
            return
         endif

         write(aval,*,iostat=ios) rval
         if(ios.ne.0) then
            call ceprint("fread_to: #Error in internal write, writing float#")
            return
         endif
         aval=adjustl(aval)
      endif
      fread_to=.true.
    end function fread_to
  end function rpf

  !**********************************************************************
!                   NoFlagInListSet
!**********************************************************************
logical function NoFlagInListSet(nev,nrfl,nlist)
  implicit none
  integer nev, nlist, nrfl(nlist)
  integer nfl
  external IsFlagSet
  logical IsFlagSet
  do nfl=1,nlist
     if(IsFlagSet(nev,nrfl(nfl))) then
        NoFlagInListSet=.false.
        return
     endif
  end do
  NoFlagInListSet=.true.
  return
end function NoFlagInListSet

!**********************************************************************
!                        get_checked_ranges
!**********************************************************************
subroutine get_checked_ranges(xu, xo, np, det, modt)
! ----------------------------------------------------------------------
! return checked axix limit 
!----------------------------------------------------------------------
  use gblcom
  implicit none
  integer, intent(in) :: np, det, modt
  real xu, xo
  real, parameter :: ulim=-1.e5, olim=1.e5  ! unteres und oberes extrem achsenlimit

  xu=xyl(np,det)
  xo=xyh(np,det)
  if(isnan(xu)) then
     xu=ulim*1.01
     call ceprint("#Previous value of lower axis limit is NAN#")
  endif
  if(isnan(xo)) then
     xo=olim*1.01
     call ceprint("#Previous value of lower axis limit is NAN#")
  endif
     

  if(xu.lt.ulim .or. xu.gt.olim) then
     print *,trim(ptx(np))//": Lower axis limit is extreme:", xu
!!$     xu=ulim
!!$     call askf("Lower axis limit: ",xu)
!!$     xyl(np,det)=xu
  endif
  if(xo.lt.ulim .or. xo.gt.olim) then
     print *,trim(ptx(np))//": Upper axis limit is extreme:",xo
!!$     xo=olim
!!$     call askf("Upper axis limit: ",xo)
!!$     xyh(np,det)=xo
  endif
  if(xu.eq.xo) then
     xo=1
     xu=-1
     xyl(np,det)=xu
     xyh(np,det)=xo
  elseif(xu.gt.xo) then
     xu=xyh(np,det)
     xo=xyl(np,det)
     xyl(np,det)=xu
     xyh(np,det)=xo
  endif

  do while(xyl(np,det).le.0. .and. modt>0)
     call ceprint("#Lower axis limit < 0. Set positive value for log scaling#")
     CALL ASKF("Lower limit:",xyl(np,det))
     xu=xyl(np,det)
     if(nast>0) return
  end do
  
  do while(xyh(np,det).le.0. .and. modt>0)
     call ceprint("#Upper limit < 0. Set positive value for log scaling#")
     CALL ASKF("Upper axis limit:",xyh(np,det))
     xo=xyh(np,det)
     if(nast>0) return
  end do
  
  return
end subroutine get_checked_ranges

end module utils

module readers
  use utils
  real, private, parameter :: dv_di=20./2**16
  integer, private :: lun, fn_was=0

contains
 
  
  subroutine rrd(data,nrec,is,filenum) 
    !************************************************************************
    !This routine is called by different threads. It remembers internally
    !which file is open from the previous call. Only one thread can access
    !at a time.
    !
    !Input:
    ! nrec = RECORDNUMMER
    ! IS   = KONTROLLVAR. : -1   = no manipulation        
    !                        0   = baseline subtraction
    !                        1   = base line subtraction + integration
    !                        2   = moving average as in cmp
    !Output:
    ! filenum  = number of file in list
    !***********************************************************************
    use gblcom
    implicit none
    integer, intent(in) :: nrec, is
    integer, intent(out) :: filenum
    real, intent(out) ::  data(krd)
    integer*2 :: sdata(krd)

    integer*8 off
    integer :: n, m, ios
    real z, yoff
    logical open_new_file

    ! find the file number
    open_new_file=.false.
    do while (nrec.gt.recsum(fn_was))
       if(nrec.gt.rectat) then
          print *,' routine rrd: requested record out of range: ',nrec
          nast=1
          return
       endif
       fn_was=fn_was+1
       open_new_file=.true.
    end do
    if(.not.open_new_file) then
       do while (nrec.le.recsum(fn_was-1))
          fn_was=fn_was-1
          open_new_file=.true.
       end do
    endif
    filenum=fn_was


    if(open_new_file) then
       close(unit=lun, iostat=ios)
       open(newunit=lun, file=trim(flist(filenum))//".rdt", &
            form="unformatted", access="stream",&
            status="old", action="read", iostat=ios)
       if(ios.ne.0) then
          call ceprint("#Routine RRD: Error opening file "//&
               trim(flist(filenum))//".rdt#")
          return 
       endif
    endif

    off=int8(nrec-1-recsum(filenum-1))*int8(krd*2 + (khd+ktd+kad)*4) + 1_8
    read(unit=lun,pos=off, iostat=ios) ih(1:khd), tt(1:ktd), dvm(1:kad), sdata
    if(ios.ne.0) then
       call ceprint("\n#Routine RRD: Error reading data: "//trim(flist(filenum))//"#")
       print *,"Record:",nrec,", offset: ",off
       print *,"First and last record in file:", recsum(filenum-1)+1, recsum(filenum)
       nast=1
       return
    end if
    ih(1)=ih(1)+1          ! detector number, DAQ counts from 0
    !******* extract digitizer data **********************************
    data=sdata*dv_di
    if(is.lt.0) return


    yoff=fip(nrec,npbase)      ! subtract baseline
    data=sdata*dv_di - yoff
    if(is.eq.0) return
    
    if(is.eq.1) then           ! and integrate
       z=0.
       do n=1,krd
          z=z+data(n)
          data(n)=z
       end do
    elseif(is.eq.2) then       ! moving average
       m=kam(ih(1))
       z=sum(data(1:m))
       do n=1,krd-m
          z=(data(n+m)-data(n))+z
          data(n)=z/m
       end do
       do n=krd-m+1,krd
          data(n)=data(krd-m)
       end do
    endif
  end subroutine rrd

  integer function read_digitizer_data(data, nrec, type)           
    !************************************************************************
    ! Works for streaming data and conventional data
    ! Does not read header data in case of conventional data
    ! This routine may be called by different threads. It remembers internally
    ! which file is open from the previous call. 
    !Input:
    ! nrec = record number
    ! type   =  -1   no manipulation        
    !           =   0   baseline subtraction
    !           =   1   baseline subtraction and integration
    !           =   2   moving average as in cmp
    ! Returns file number from which data were read or -1 in case of error
    ! In case of error nast is set to 1
    !***********************************************************************
    use gblcom, only : nast, krd, kam, npbase, npdet,&
          use_ofstamps, ip, fip
    implicit none
    integer, intent(in) :: nrec, type 
    real, intent(out) ::  data(krd) 

    integer :: n, m
    real :: z, rinvm

    if(use_ofstamps) then
       read_digitizer_data=read_streaming_data(nrec,data)
    else
       read_digitizer_data=read_conventional_data(nrec,data)
    endif

    if(read_digitizer_data.lt.0) then
       nast=1
       return
    endif


    if(type.lt.0) return            ! raw

    z=fip(nrec,npbase)             
    data(:) = data(:) - z
    if(type.eq.0) return            ! subtract base line


    if(type.eq.1) then              ! integration
       z=0.
       do n=1,krd
          z=z+data(n)
          data(n)=z
       end do
    else if(type.eq.2) then              ! moving average
       m=kam(ip(nrec,npdet))
       rinvm=1./float(m)
       z=sum(data(1:m))
       do n=1,krd-m
          z=(data(n+m)-data(n))+z
          data(n)=z*rinvm
       end do
       z=data(krd-m)
       data(krd-m+1:krd)=z
    endif
  end function read_digitizer_data


  !**********************************************************************
  !                      read_conventioal_data
  !**********************************************************************
  integer function read_conventional_data(rec, data)
    use gblcom, only : krd, khd, ktd, kad, npfile, flist, recsum, ip
    implicit none
    integer, intent(in) :: rec
    real, intent(out) :: data(krd)
    integer*8 :: bytes_per_header, bytes_per_record, off
    integer filenum, ios
    integer*2 :: sdata(krd)

    read_conventional_data=-1
    filenum=ip(rec,npfile)
    if(fn_was.ne.filenum) then
       close(unit=lun, iostat=ios)
       open(newunit=lun, file=trim(flist(filenum))//".rdt", form="unformatted", &
            access="stream", status="old",  position="rewind", &
            action="read", iostat=ios)
       if(ios.ne.0) then
          call ceprint("\n#Error opening file in routine read_conventional data#")
          call ceprint("#file: "//trim(flist(filenum))//".rdt#")
          return 
       endif
       fn_was=filenum
    endif
    read_conventional_data=filenum        
    bytes_per_header = (khd+ktd+kad)*4
    bytes_per_record = krd*2 + bytes_per_header

    off = 1_8 + int8(rec-1-recsum(filenum-1))*bytes_per_record + bytes_per_header
    read(unit=lun,pos=off, iostat=ios) sdata
    if(ios.ne.0) then
       call ceprint("\n#Error in routine read_conventional_data#")
       call ceprint("#file: "//trim(flist(filenum))//"#")
       call cenprint("#Record in file: #")
       print*,rec-recsum(filenum-1)
       return
    endif
    data=sdata*dv_di
    return
  end function read_conventional_data


  !**********************************************************************
  !                 read_streaming_data
  !**********************************************************************
  integer function read_streaming_data(rec, data)
    ! returns filenumber or -1 in case of error
    use stream_helpers, only : open_stream_file,&
         close_this_stream, get_stream_file, get_lun_stream, get_stream_channel
    use gblcom, only : krd, ndetm, npfile, npdet, pre_trigger, time_base,&
         stamp, ip, flist, stream_delay
    use ISO_C_BINDING, only : C_INT64_T
    implicit none
    integer, intent(in) :: rec
    real, intent(out) :: data(krd)

    integer*8 :: off, time_base_i8, npre
    integer*8,save :: file_size(ndetm)
    integer, save :: filenum_was(ndetm)
    integer det, filenum, stream_channel, ios
    integer*2 sdata(krd)
    character(len=80) temp_str
    
    external get_stream_delay   ! c-function get_stream_delay_
    integer(kind=C_INT64_T) get_stream_delay
    
    read_streaming_data=-1
    filenum=ip(rec,npfile)
    det=ip(rec,npdet)
    stream_channel=get_stream_channel(det)
    if(stream_channel.lt.0) then
       call ceprint("#Fatal error: mapping of detectors to streaming channels not defined#")
       stop
    endif
    ! isticks i1,i2 and i3 access same stream file. Thats why file size
    ! is labeled with stream_channel instead of detector number 

  
    if(filenum .ne. filenum_was(stream_channel)) then
       call close_this_stream(det, filenum_was(stream_channel))
       file_size(stream_channel)=open_stream_file(det,filenum,ios)
       if(ios .ne. 0) then
          call ceprint("\n #Error opening "//trim(get_stream_file(det,filenum))//&
               " in read_streaming_data#")
          print "(2(A,i4))", " detector: ", det, ", stream channel: ", stream_channel
          print "(A,i4)"," file number:",filenum
          stop
       endif
      
       stream_delay(filenum)=get_stream_delay(&
            trim(flist(filenum))//".dig_stamps"//char(0))
       filenum_was(stream_channel)=filenum
    endif

    time_base_i8=10*nint(time_base)
    npre=pre_trigger*krd/8

    off=((stamp(rec) - stream_delay(filenum))/time_base_i8 - npre)*2_8 +1_8
    off=max(off, 1_8)
    if(off+2*krd .gt. file_size(stream_channel)) then
       call ceprint("\n# Error in routine read_streaming_data#")
       call ceprint("#file:"//trim(get_stream_file(det,filenum))//"#\n")
       call ceprint("#Requested to read beyond file end#")
       call cenprint("#Offset shifted left [samples]:#")
       return
    endif

    read(unit=get_lun_stream(det), pos=off, iostat=ios) sdata
    if(ios .ne. 0) then
       call cenprint("\n# Error in routine read_streaming_data, record:#")
       print *,rec
       call ceprint("#File:"//trim(get_stream_file(det,filenum))//"#")
       write(temp_str,"(3(A,i4))") "file number: ", filenum," detector: ", &
            det, ", stream channel: ", stream_channel
       call ceprint("#"//trim(temp_str)//"#")
       print *,"(file_size-offset)/2", (file_size(stream_channel)-off)/2
       print *,"offset",off
       return
    endif
    data(:)=sdata(:)*dv_di
    read_streaming_data=filenum
  end function read_streaming_data

end module readers

