# This script is meant to be sourced in a test script to start a docker container to test in
# it interprets the options:
# -i|-image|--image: docker image to start from
# -a|-arch|--arch: 64, x86_64 or linux-x86_64 for 64 bit Linux build (default); ix86, 32 or linux-ix86 for 32 bits Linux build; win, windows-x86_64 or mingw-w64 for Windows 64 bit build
# -b|-bits|--bits: select 32 or 64 bits Linux build (default 64 = same as -arch x86_64)
# -d|-bindir|--bindir: test directory (default ~/bin)
# -d|-testdir|--testdir: test directory (default ~/genomecomb.smalltestdata)
# and presents all options (excluding -bbuilddir) to the program in the docker

# you can use the following to use in the test script
## Prepare and start docker
## ========================
#script="$(readlink -f "$0")"
#dir="$(dirname "$script")"
#source "${dir}/hbb_start.sh"

# This test environment requires docker, make sure it is installed
# e.g. on ubuntu and derivatives
# sudo apt install docker.io
# Also make sure you have permission to use docker
# sudo usermod -a -G docker $USER

# Prepare and start docker
# ========================

# check if we are already in test box
if [ ! -f /io/tests/dockertest.sh ]; then
	# find directory of script
	script="$(readlink -f "$0")"
	dir="$(dirname "$script")"
	file="$(basename "$script")"
	if [ $(basename "$dir") = "tests" ] ; then
		srcdir="$(dirname "$dir")"
		file="tests/$file"
	else
		srcdir="$dir"
	fi

	echo "Running script $dir/$file"
	testdir=""
	bindir=""
	image=centos:7
	arch=linux-x86_64
	argumentspos=1; 
	while [[ "$#" -gt 0 ]]; do case $1 in
		-i|-image|--image) image=$2; shift;;
		-b|-bits|--bits) arch="$2"; shift;;
		-a|-arch|--arch) arch="$2"; shift;;
		-d|-testdir|--testdir) testdir="$(readlink -f "$2")" ; shift;;
		-d|-bindir|--bindir) bindir="$(readlink -f "$2")" ; shift;;
		*) arguments[$argumentspos]="$1"; argumentspos+=1 ; arguments[$argumentspos]="$2"; argumentspos+=1 ; shift;;
	esac; shift; done
	if [ "$bits" = "32" ] ; then arch=linux-ix86 ; fi
	if [ "$bits" = "64" ] ; then arch=linux-x86_64 ; fi
	if [ "$arch" = "32" ] || [ "$arch" = "ix86" ] ; then 
		arch=linux-ix86
	elif [ "$arch" = "64" ] || [ "$arch" = "x86_64" ]; then
		arch=linux-x86_64
	elif [ "$arch" = "win" ] || [ "$arch" = "mingw-w64" ]; then
		arch=windows-x86_64
	fi
	if [ "$testdir" = "" ] ; then
		testdir="$HOME/genomecomb.smalltestdata"
	fi
	if [ "$bindir" = "" ] ; then
		bindir="$HOME/bin"
	fi
	echo "test in $image"
	echo "testdir=$testdir"
	echo "bindir=$bindir"
	echo "srcdir=$srcdir"
	# run the script in test container
	uid=$(id -u)
	gid=$(id -g $uid)
	
	if [ "$arch" = "linux-ix86" ] ; 	then
		docker run --net=host -t -i --rm -v "$srcdir:/io" -v "$testdir:/test" -v "$bindir:/testbin" "$image" linux32 bash "/io/$file" "stage2" "$file" "$arch" "$uid" "$gid" "$srcdir" "$testdir" "$bindir" "$image" ${arguments[*]}
	else
		docker run --net=host -t -i --rm -v "$srcdir:/io" -v "$testdir:/test" -v "$bindir:/testbin" "$image" bash "/io/$file" "stage2" "$file" "$arch" "$uid" "$gid" "$srcdir" "$testdir" "$bindir" "$image" ${arguments[*]}
	fi
	exit
fi

if [ "$1" = "stage2" ] ; then
	# in stage 2 we will create the user test with sudo access
	# then restart the script (skipping to stage 3)
	file=$2
	arch=$3
	uid=$4
	gid=$5
	srcdir=$6
	testdir=$7
	bindir=$8
	image=$9
	echo "preparing user test with uid=$uid and gid=$gid"
	groupadd test --gid $gid
	useradd test -m --uid $uid --gid $gid
	# prepare the user test with sudo rights
	if [[ "$image" =~ centos ]] ; then
		echo "installing sudo ($image)"
		# to stop "checksum is invalid" errors when using yum in 32 bit docker
		if [ "$image" =~ linux-ix86 ] ; then
			if ! rpm --quiet --query yum-plugin-ovl; then
				yum install -q -y yum-plugin-ovl
			fi
		fi
		if ! rpm --quiet --query sudo; then
			yum install -q -y sudo
		fi
		# usermod -a -G wheel test
		echo "test ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-test
	elif [[ "$image" =~ ubuntu ]] ; then
		echo "installing sudo ($image)"
		apt update
		apt -y install sudo
		adduser test sudo
		echo "test ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-test
	else
		echo "did not install sudo ($image)"
	fi
	# (re)start script for stage 3: running the actual code
	sudo -u test bash /io/$file "stage3" "$file" "$bits" "$uid" "$gid" "$srcdir" "$testdir" "$bindir" "$image" ${arguments[*]}
	exit
fi

# stage 3: run the actual script (first do some settings)

function yuminstall {
	echo "yuminstall $1"
	if ! rpm --quiet --query "$1"; then
		sudo yum install -y "$1"
	fi
}

file=$2

if [ $(basename "$file") = "dockertest.sh" ] ; then
	# install yuminstall in .bashrc so it will be available in the new shell started here
	# sudo mkdir -p /home/test
	cd /home/test
	ln -s /test /home/test/genomecomb.smalltestdata
	ln -s /testbin /home/test/bin
	echo "
	file=$2
	arch=$3
	uid=$4
	gid=$5
	srcdir=$6
	testdir=$7
	bindir=$8
	image=$9
	" >> /home/test/.bashrc
	echo 'if [ "$arch" = 'linux-ix86' ] ; then
		ARCH='-linux-ix86'
		arch=linux-ix86
		bits=32
	else
		ARCH=''
		arch=linux-x86_64
		bits=64
	fi
	function yuminstall {
		echo "yuminstall $1"
		if ! rpm --quiet --query "$1"; then
			sudo yum install -y "$1"
		fi
	}
	export PATH=/home/test/bin:$PATH
	' >> /home/test/.bashrc
	# if run as dockertest.sh directly, show a shell
	cd /home/test
	echo "shell started by dockertest.sh"
	bash
	exit
fi

# if sourced in another script, continue executing this other script
if [ "$3" = 'linux-ix86' ] ; then
	ARCH='-linux-ix86'
	arch=linux-ix86
	bits=32
else
	ARCH=''
	arch=linux-x86_64
	bits=64
fi
uid=$4;
gid=$5;
srcdir=$6;
testdir=$7;
bindir=$8;
image=$9;
shift 9;

echo "Entering docker environment; testing using image $image, arch $arch, $bits bits, uid=$uid, gid=$gid"
