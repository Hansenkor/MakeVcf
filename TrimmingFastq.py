"""
    date: 2025.02.28
    object: Trimming.py (FASTQ 파일들에 대해 Trimmomatic trimming 분석을 자동화하는 스크립트)
    author: hansen

This script requires that Java is installed and that Trimmomatic is accessible via the provided jar file.
Usage:
  TrimmingFastq.py [--work_dir=<dir>] [--in_dir=<dir>] [--out_dir=<dir>] [--jar_path=<jar>] [--threads=<n>] [--adapter_fa=<path>] [--smm=<val>] [--pal=<val>] [--thr=<val>] [--samples=<list>] [--lead=<val>] [--trail=<val>] [--slidewin=<val>] [--minlen=<val>]

Options:
  --work_dir=<dir>      작업 디렉토리.
  --in_dir=<dir>        입력 FASTQ 파일들이 위치한 디렉토리.
  --out_dir=<dir>       결과가 저장될 출력 디렉토리. 지정되지 않으면 workdif/trimmed 생성됨.
  --jar_path=<jar>      Trimmomatic jar 파일의 경로. [default: /data/share/Trimmomatic-0.39/trimmomatic-0.39.jar]
  --threads=<n>         사용할 스레드 수. [default: 30]
  --adapter_fa=<path>   ILLUMINACLIP에 사용할 adapter 파일에 얼마나 서열을  [default: /data/share/Trimmomatic-0.39/adapters/TruSeq3-PE.fa]
  --smm=<val>           seed mismatches로 허용가능한 seed 서열의 mismatch 개수. [default: 2]
  --pal=<val>           palindrome clip threshold로 paired end에 해당하는 경우 target read 앞뒤에 존재하는 adapter 서열을 역상보적으로 삭제 [default: 30] 
  --thr=<val>           simple clip threshold로 설정기준에 부합하는 충분히 정확한 match가 확인되면 적절히 clipping함. [default: 10] 
  --samples=<list>      샘플 리스트 (따옴표(" ") 내에 공백으로 구분된 샘플 이름들). ex: "A-1 A-2 B-1 B-2 B-3 B-4 B-5 C-1 C-2"
  --lead=<val>          LEADING 옵션 값. [default: 3]
  --trail=<val>         TRAILING 옵션 값. [default: 3]
  --slidewin=<val>      SLIDINGWINDOW 옵션 값 (window:quality). [default: 4:15]
  --minlen=<val>        MINLEN 옵션 값. [default: 36]

  <sample>_1.fastq.gz <sample>_2.fastq.gz [ex : A-1_1.fastq.gz , A-1_2.fastq.gz ] 이런식으로 파일명이 되어있어야 --samples의 샘플리스트를 읽을 수 있습니다.
"""


import os
import sys
import logging
import subprocess
from pytz import timezone
from docopt import docopt# conda install -c conda-forge docopt 해야 
from datetime import datetime



def make_logger(log_dir, name="app", consoleset=True, streamset=True) -> logging.Logger:
    """
    Logger 생성 함수.
    
    - 콘솔 핸들러와 파일 핸들러를 추가하여, 로그 메시지를 콘솔 및 파일에 기록합니다.
    - 로그 시간은 Asia/Seoul 시간대를 사용합니다.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    loggerformat = "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(loggerformat)
        
    #  change time zone 시간대 설정
    tz = timezone("Asia/Seoul")
    def timetz(*args):
        return datetime.now(tz).timetuple()
    logging.Formatter.converter = timetz

    if consoleset:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(formatter)
        logger.addHandler(console)
    if streamset:
        loggerfile = os.path.join(log_dir, name)
        file_handler = logging.FileHandler(filename=f"{loggerfile}.log")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger



class TrimmingFastq:
    @staticmethod
    def check_java(logger):
        """Java가 설치되어 있는지 확인합니다."""
        try:
            subprocess.run(["java", "-version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            logger.info("Java 설치 확인 완료.")
        except subprocess.CalledProcessError:
            logger.error("Java가 설치되어 있지 않거나 실행할 수 없습니다.")
            sys.exit(1)

    @staticmethod
    def check_jar(jar_path: str, logger) -> None:
        """Trimmomatic jar 파일이 존재하는지 확인합니다."""
        if not os.path.isfile(jar_path):
            logger.error(f"Jar 파일이 존재하지 않습니다: {jar_path}\nJar 파일의 경로를 확인하거나, 올바른 경로를 지정해 주세요.")
            sys.exit(1)
        logger.info(f"Jar 파일 확인: {jar_path}")

    
    @staticmethod
    def check_in_dir(in_dir: str, logger) -> None:
        """입력 디렉토리가 존재하는지 확인합니다."""
        if not os.path.isdir(in_dir):
            logger.error(f"입력 디렉토리가 존재하지 않습니다: {in_dir}\nin_dir 경로를 확인해 주세요.")
            sys.exit(1)
        logger.info(f"입력 디렉토리 확인: {in_dir}")

    @staticmethod
    def check_workdir(work_dir: str, logger) -> None:
        """작업 디렉토리가 존재하는지 확인합니다."""
        if not os.path.isdir(work_dir):
            TrimmingFastq.logger.error(f"작업 디렉토리가 존재하지 않습니다: {work_dir}\nwork_dir 경로를 확인해 주세요.")
            sys.exit(1)
        logger.info(f"작업 디렉토리 확인: {work_dir}")

    @staticmethod
    def check_sample_files(in_dir: str, samples_str: str, logger) -> None:
        """샘플 목록에 따른 FASTQ 파일들이 존재하는지 미리 확인합니다."""
        missing_samples = []
        for sample in samples_str.split():
            in_file1 = os.path.join(in_dir, f"{sample}_1.fastq.gz")
            if not os.path.isfile(in_file1):
                missing_samples.append(sample)
        if missing_samples:
            logger.error(f"다음 샘플의 FASTQ 파일이 존재하지 않습니다: {', '.join(missing_samples)}")
            sys.exit(1)
        logger.info("모든 지정된 샘플의 FASTQ 파일 확인 완료.")


    def get_samples(in_dir: str, samples_arg: str = None) -> list:
        """
        --samples 옵션이 제공된 경우, 해당 문자열을 공백 기준으로 분할하여 리스트 반환.
        제공되지 않은 경우, in_dir 내의 FASTQ 파일명을 통해 샘플 리스트를 생성합니다.
        """
        if samples_arg:
            TrimmingFastq.check_sample_files(in_dir, samples_arg)
            samples = samples_arg.split()
            logger.info(f"사용자 지정 샘플 리스트: {samples}")
            return samples
        else:
            sample_files = [f for f in os.listdir(in_dir)
                            if f.endswith("_1.fastq.gz") or f.endswith("_1_paired.fastq.gz")]
            if not sample_files:
                logger.error("FASTQ 파일의 형태가 sample 표시에 적합하지 않습니다.")
                sys.exit(1)
            
            samples = set()
            for f in sample_files:
                if f.endswith("_1_paired.fastq.gz"):
                    sample = f.replace("_1_paired.fastq.gz", "")
                else:
                    sample = f.replace("_1.fastq.gz", "")
                samples.add(sample)
            samples_list = list(samples)
            logger.info(f"자동 생성된 샘플 리스트: {samples_list}")
            return samples_list
    
    
    @staticmethod
    def check_adapter(adapter_fa: str, logger) -> None:
        """Adapter 파일이 존재하는지 확인합니다."""
        if not os.path.isfile(adapter_fa):
            TrimmingFastq.logger.error(f"Adapter 파일이 존재하지 않습니다: {adapter_fa}\nAdapter 파일 경로를 확인하거나, 올바른 경로를 지정해 주세요.")
            sys.exit(1)
        logger.info(f"Adapter 파일 확인: {adapter_fa}")

    


    @staticmethod
    def run_trimmomatic(work_dir: str, in_dir: str, out_dir: str, jar_path: str, threads: int,adapter_fa: str, smm: int, pal: int, thr: int, samples: list,lead: int, trail: int, slidewin: str, minlen: int, logger) -> None:
        """
        각 샘플에 대해 Trimmomatic을 실행합니다.
          - PE: {in_dir}/{SAMPLE}_1.fastq.gz와 {in_dir}/{SAMPLE}_2.fastq.gz가 모두 존재하면,
                {out_dir}/{SAMPLE}_1_paired.fastq.gz, {out_dir}/{SAMPLE}_1_unpaired.fastq.gz,
                {out_dir}/{SAMPLE}_2_paired.fastq.gz, {out_dir}/{SAMPLE}_2_unpaired.fastq.gz 생성.
          - SE: {in_dir}/{SAMPLE}_1.fastq.gz만 존재하면,
                {out_dir}/{SAMPLE}_trimmed.fastq.gz 생성.
        """
        if not out_dir:
            out_dir = os.path.join(work_dir, "trimmed")
        os.makedirs(out_dir, exist_ok=True)
        logger.info(f"출력 디렉토리 생성/확인: {out_dir}")

        try:
            os.chdir(work_dir)
            logger.info(f"작업 디렉토리 변경: {work_dir}")
        except Exception as exc:
            logger.error(f"작업 디렉토리 변경 중 에러 발생: {exc}")
            sys.exit(1)

        for sample in samples:
            in_file1 = os.path.join(in_dir, f"{sample}_1.fastq.gz")
            in_file2 = os.path.join(in_dir, f"{sample}_2.fastq.gz")
            
            # Single-End (SE) 처리: _1.fastq.gz만 있을 경우
            if os.path.isfile(in_file1) and not os.path.isfile(in_file2):
                out_file = os.path.join(out_dir, f"{sample}_trimmed.fastq.gz")
                command = (
                    f"java -jar {jar_path} SE -threads {threads} -phred33 "
                    f"{in_file1} {out_file} "
                    f"ILLUMINACLIP:{adapter_fa}:{smm}:0:{thr} LEADING:{lead} TRAILING:{trail} "
                    f"SLIDINGWINDOW:{slidewin} MINLEN:{minlen}"
                )
                TrimmingFastq.logger.info(f"Trimmomatic (SE) 실행 중: {sample}")
                try:
                    subprocess.run(command, shell=True, check=True)
                    TrimmingFastq.logger.info(f"{sample} (SE) 처리 완료.")
                except subprocess.CalledProcessError as error:
                    TrimmingFastq.logger.error(f"Trimmomatic 실행 중 에러 발생 ({sample}, SE): {error}")
                    continue

            # Paired-End (PE) 처리: _1.fastq.gz와 _2.fastq.gz 모두 존재하는 경우
            elif os.path.isfile(in_file1) and os.path.isfile(in_file2):
                out1_paired = os.path.join(out_dir, f"{sample}_1_paired.fastq.gz")
                out1_unpaired = os.path.join(out_dir, f"{sample}_1_unpaired.fastq.gz")
                out2_paired = os.path.join(out_dir, f"{sample}_2_paired.fastq.gz")
                out2_unpaired = os.path.join(out_dir, f"{sample}_2_unpaired.fastq.gz")
                
                command = (
                    f"java -jar {jar_path} PE -threads {threads} -phred33 "
                    f"{in_file1} {in_file2} "
                    f"{out1_paired} {out1_unpaired} {out2_paired} {out2_unpaired} "
                    f"ILLUMINACLIP:{adapter_fa}:{smm}:{pal}:{thr} LEADING:{lead} TRAILING:{trail} "
                    f"SLIDINGWINDOW:{slidewin} MINLEN:{minlen}"
                )
                TrimmingFastq.logger.info(f"Trimmomatic (PE) 실행 중: {sample}")
                try:
                    subprocess.run(command, shell=True, check=True)
                    TrimmingFastq.logger.info(f"{sample} (PE) 처리 완료.")
                except subprocess.CalledProcessError as error:
                    TrimmingFastq.logger.error(f"Trimmomatic 실행 중 에러 발생 ({sample}, PE): {error}")
                    continue

            else:
                TrimmingFastq.logger.error(f"{sample}의 입력 파일이 올바르지 않습니다: {in_file1} 혹은 {in_file2}")
                continue

        TrimmingFastq.logger.info("All done!")


    """
    --SE                  Single-end read input. Default input choice is single-
                        end if nothing is specified
    --PE                  Paired-end read input. Must have the exact same file
                        name and end with _F for the forward read and _R for
                        the reverse read


    Trimmomatic options:
    The order the options will be run are: ILLUMINACLIP, LEADING,
    TRAILING, SLIDINGWINDOW, MINLEN

    --clip=ILLUMINACLIP
                        ILLUMINACLIP options. MiSeq & HiSeq usually
                        TruSeq3.fa; GAII usually TruSeq2.fa. Default is
                        ILLUMINACLIP:TruSeq3-SE.fa:2:30:10. Usage:
                        --clip=<adapterseqs>:<seed mismatches>:<palindrome
                        clip threshold>:<simple clip threshold>
    --lead=LEADING      Set the minimun quality required to keep a base.
                        Default is LEADING=3. Usage: --lead=<quality>
    --trail=TRAILING    Set the minimum quality required to keep a base.
                        Default is TRAILING=3. Usage: --trail=<quality>
    --slidewin=SLIDINGWINDOW
                        SLIDINGWINDOW options. Default is SLIDINGWINDOW:4:15.
                        Usage: --slidewin=<window_size>:<required_quality>

    HTSeq options:
    --stranded=STRANDED
                        Stranded options: yes, no, reverse. Default is
                        --stranded=reverse. Usage: --stranded=yes/no/reverse
    --order=ORDER       Order options: name, pos. Usage: --order=name/pos.
    --minqual=MINQUAL   Skip all reads with quality lower than the given
                        value. Default is --minqual=10. Usage:
                        --minqual=<value>
    --idattr=IDATTR     Feature ID from the GTF file to identify counts in the
                        output table Default is --idattr=gene_id. Usage:
                        --idattr=<id attribute>
    --mode=MODE         Mode to handle reads overlapping more than one
                        feature. Default is --mode=union. Usage: --mode=union
                        /intersection-strict/intersection-nonempty

    """
    @staticmethod
    def main() -> None:
        arguments = docopt(__doc__)
        # work_dir가 제공되지 않으면 현재 디렉토리를 기본값으로 사용합니다.
        work_dir = arguments.get("--work_dir")
        if not work_dir:
            work_dir = os.getcwd()
        # --in_dir 옵션이 지정되지 않으면 원본 FASTQ 디렉토리 (예: work_dir/fastq) 사용
        in_dir = arguments.get("--in_dir")
        if not in_dir:
            in_dir = os.path.join(work_dir, "fastq")
            
        out_dir = arguments["--out_dir"]
        if not out_dir: 
            out_dir = os.path.join(work_dir,"trimmed")
        # work_dir = arguments.get("--work_dir") or os.getcwd()
        # in_dir = arguments.get("--in_dir") or os.path.join(work_dir, "fastq")
        # out_dir = arguments.get("--out_dir") or os.path.join(work_dir, "trimmed")
        jar_path = arguments["--jar_path"]
        threads = int(arguments["--threads"])
        adapter_fa = arguments["--adapter_fa"]
        smm = int(arguments["--smm"])
        pal = int(arguments["--pal"])
        thr = int(arguments["--thr"])
        samples_arg = arguments.get("--samples")
        samples = arguments["--samples"]
        lead = int(arguments["--lead"])
        trail = int(arguments["--trail"])
        slidewin = arguments["--slidewin"]
        minlen = int(arguments["--minlen"])
        

        # logger 생성 및 클래스 변수에 할당 (로그는 work_dir에 저장됨)
        logger = make_logger(work_dir, name="TrimmingFastq")
        TrimmingFastq.logger = logger
        TrimmingFastq.logger.info("TrimmingFastq 분석 시작.")

        # 필수 도구 및 파일 존재 여부 확인
        TrimmingFastq.check_java()
        TrimmingFastq.check_adapter(adapter_fa)
        TrimmingFastq.check_jar(jar_path)
        TrimmingFastq.check_workdir(work_dir)
        TrimmingFastq.check_in_dir(in_dir)
                
        # get_samples를 호출하여 샘플 리스트를 얻음 (문자열 옵션이든 자동 생성이든)
        samples = TrimmingFastq.get_samples(in_dir, samples_arg)

        TrimmingFastq.run_trimmomatic(work_dir, in_dir, out_dir, jar_path, threads, adapter_fa, smm, pal, thr, samples, lead, trail, slidewin, minlen)


if __name__ == "__main__":
    TrimmingFastq.main()