#!/usr/bin/env python3
"""Preflight for the narrowly supported floating-point reproduction target.

This probe does not execute the global certificate or establish its numerical
premises. It deliberately rejects wider long-double formats: the frozen
extension's 30-digit decimal writers were reviewed for a 64-bit significand.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import tempfile

FLAGS = ["-std=c++17", "-O2", "-DNDEBUG", "-fno-fast-math",
         "-ffp-contract=off", "-frounding-math"]
PROBE = r'''
#include <cfloat>
#include <cfenv>
#include <cstdint>
#include <iostream>
#include <limits>
int main() {
  const std::uint16_t one=1;
  std::cout << "{\"radix\":" << FLT_RADIX
    << ",\"double_bytes\":" << sizeof(double)
    << ",\"double_significand_bits\":" << DBL_MANT_DIG
    << ",\"double_max_exponent\":" << DBL_MAX_EXP
    << ",\"double_iec559\":" << std::numeric_limits<double>::is_iec559
    << ",\"long_double_significand_bits\":" << LDBL_MANT_DIG
    << ",\"long_double_max_exponent\":" << LDBL_MAX_EXP
    << ",\"long_double_bytes\":" << sizeof(long double)
    << ",\"round_to_nearest\":" << (std::fegetround()==FE_TONEAREST)
    << ",\"little_endian\":"
    << (*reinterpret_cast<const unsigned char*>(&one)==1)
    << "}\n";
}
'''


def validate_environment(info):
    expected = {
        "radix": 2, "double_bytes": 8, "double_significand_bits": 53,
        "double_max_exponent": 1024, "double_iec559": 1,
        "long_double_significand_bits": 64, "long_double_max_exponent": 16384,
        "round_to_nearest": 1, "little_endian": 1,
    }
    for key, value in expected.items():
        if info.get(key) != value:
            raise ValueError(f"Unsupported reproduction environment: {key}="
                             f"{info.get(key)!r}; expected {value}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="xy-environment-") as tmp:
        root = Path(tmp)
        source, binary = root/"probe.cpp", root/"probe"
        source.write_text(PROBE)
        subprocess.run([args.compiler, *FLAGS, str(source), "-o", str(binary)],
                       check=True, capture_output=True, text=True)
        output = subprocess.run([str(binary)], check=True, capture_output=True,
                                text=True)
        info = json.loads(output.stdout)
        validate_environment(info)
        compiler = subprocess.run([args.compiler, "--version"], check=True,
                                  capture_output=True, text=True).stdout.splitlines()[0]
    result = {
        "status": "SUPPORTED_EXTENDED64_ENVIRONMENT_PASS",
        "scope": "Environment preflight only; not the global certificate",
        "global_certificate_reexecuted": False,
        "compiler": compiler, "flags": FLAGS,
        "platform": platform.platform(), "probe": info,
        "probe_source_sha256": hashlib.sha256(PROBE.encode()).hexdigest(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
