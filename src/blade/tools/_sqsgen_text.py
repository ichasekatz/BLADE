"""Build the ``sqsgen.in`` composition text for ATAT mcsqs runs.

Standalone function extracted from :class:`~blade.tools.blade_sqsgen.BladeSQS`
that constructs the ``sqsgen.in`` file content from a list of composition level
seeds.  Imported by ``blade_sqsgen`` — do not call directly unless you know
what you are doing.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import gcd

__author__ = "Chase Katz"

# sqsgen.in level indices used for binary systems (len_comp == 2).
# Levels 0, 1, 2, and 5 cover endmembers, 50/50, and the primary binary
# compositions in the standard BLADE sqsgen_levels seed list.
_binary_level_indices: tuple[int, ...] = (0, 1, 2, 5)


def _build_sqsgen_text(
    sqsgen_levels: list[dict],
    level: int,
    len_comp: int,
    unit_cell: str,
    fixed_compositions: dict[str, list[float]],
) -> str:
    """Build the ``sqsgen.in`` content string from *sqsgen_levels*.

    Iterates over the applicable level indices (determined by *level* and
    *len_comp*), deduplicates compositions by keeping the lowest level index
    for each unique composition tuple, and emits one
    ``level=N  <letter>=<comp> ...`` line per entry.

    Multi-sublattice lines are generated as the Cartesian product of all
    non-endmember compositions across active sublattice letters.
    Fixed-composition sublattices (from *fixed_compositions*) are appended
    verbatim to every emitted line.

    Args:
        sqsgen_levels (list[dict]): Ordered list of composition level
            definitions for ``sqsgen.in``.  Each entry must contain:

            - ``"level"`` (int): Level index.
            - ``"letter"`` (list[str]): Sublattice letter(s).
            - ``"compositions"`` (list[list[float]]): Fractional compositions
              for that sublattice.

        level (int): Maximum sqsgen level to include (inclusive).
        len_comp (int): Number of elements in the chemical system.
        unit_cell (str): Fractional coordinate block string (one atom per
            line) used to discover which sublattice letters are variable.
        fixed_compositions (dict[str, list[float]]): Maps sublattice letters
            to their fixed (constant) compositions.  These are appended to
            every emitted line rather than being included in the Cartesian
            product expansion.

    Returns:
        str: Complete ``sqsgen.in`` file content, or an empty string if
        no variable sublattices exist and no fixed compositions are set.
    """

    def _trim_comp(vals: list[float]) -> list[float]:
        """Pad *vals* to ``len_comp`` with zeros, then strip trailing zeros.

        Args:
            vals (list[float]): Raw composition fractions.

        Returns:
            list[float]: Trimmed composition list (at least length 1).
        """
        vals = list(vals)
        vals = vals + [0.0] * max(0, len_comp - len(vals))
        while len(vals) > 1 and vals[-1] == 0.0:
            vals.pop()
        return vals

    def _fmt_comp(vals: list[float]) -> str:
        """Format a composition list as a comma-separated string.

        Args:
            vals (list[float]): Composition fractions.

        Returns:
            str: Comma-joined string, e.g. ``"0.5,0.5"``.
        """
        return ",".join(map(str, vals))

    def _is_pure_one(vals: list[float]) -> bool:
        """Return ``True`` if *vals* represents a pure endmember (``[1.0]``).

        Args:
            vals (list[float]): Composition fractions.

        Returns:
            bool: ``True`` when *vals* equals ``[1.0]``.
        """
        return vals == [1.0]

    def _lcm(a: int, b: int) -> int:
        """Return the least common multiple of *a* and *b*.

        Args:
            a (int): First integer.
            b (int): Second integer.

        Returns:
            int: Least common multiple.
        """
        return a * b // gcd(a, b)

    def _composition_denominator(vals: list[float]) -> int:
        """Return the LCM of all fraction denominators in *vals*.

        Uses :class:`fractions.Fraction` with a limit denominator of 100
        to convert floats to exact rationals before computing the LCM.

        Args:
            vals (list[float]): Composition fractions.

        Returns:
            int: Shared denominator (e.g., ``4`` for ``[0.75, 0.25]``).
        """
        denom = 1
        for f in vals:
            frac = Fraction(float(f)).limit_denominator(100)
            denom = _lcm(denom, frac.denominator)
        return denom

    var_letters: set[str] = set()

    for line in unit_cell.strip().splitlines():
        parts = line.strip().split()
        if len(parts) >= 4:
            label = parts[3]
            if label.islower() and len(label) == 1:
                var_letters.add(label)

    sqsgen = ""

    if level >= 1 and len_comp == 1:
        indices = [0]
    elif level >= 3 and len_comp == 2:
        indices = list(_binary_level_indices)
    else:
        indices = list(range(level + 1))

    best_level_for_comp: dict[tuple[float, ...], int] = {}

    for i in indices:
        level_def = sqsgen_levels[i]
        level_num = level_def["level"]

        for comp in level_def["compositions"]:
            non_zero = [float(f) for f in comp if float(f) > 0.0]

            if non_zero == [1.0]:
                vals = _trim_comp(list(comp))
                key = tuple(vals)

                if key not in best_level_for_comp or level_num < best_level_for_comp[key]:
                    best_level_for_comp[key] = level_num

                continue

            vals = _trim_comp([float(f) for f in comp])
            key = tuple(vals)

            if key not in best_level_for_comp or level_num < best_level_for_comp[key]:
                best_level_for_comp[key] = level_num

    comp_entries: list[tuple[int, list[float]]] = [(level_num, list(comp)) for comp, level_num in best_level_for_comp.items()]

    comp_entries = sorted(comp_entries, key=lambda x: (x[0], x[1]))

    all_letters: list[str] = []

    for i in indices:
        level_def = sqsgen_levels[i]
        for letter in level_def["letter"]:
            # Skip letters with a fixed composition — they are appended separately
            if letter in fixed_compositions:
                continue
            if letter in var_letters and letter not in all_letters:
                all_letters.append(letter)

    if not all_letters:
        # All variable sublattices are fixed — emit a single level=0 line
        # so sqs2tdb -mk can create the sqsdb directory.
        if fixed_compositions:
            fixed_suffix = "".join(f"\t\t{letter}={_fmt_comp(comp)}" for letter, comp in fixed_compositions.items())
            return f"level=0{fixed_suffix}\n"
        return ""

    if len(all_letters) == 1:
        letter = all_letters[0]

        for level_num, vals in comp_entries:
            if level_num > level:
                continue

            line = f"level={level_num}\t\t{letter}={_fmt_comp(vals)}"
            if fixed_compositions:
                line += "".join(f"\t\t{fl}={_fmt_comp(fc)}" for fl, fc in fixed_compositions.items())
            sqsgen += line + "\n"

        return sqsgen

    endmember_line = "level=0"
    for letter in all_letters:
        endmember_line += f"\t\t{letter}=1.0"
    sqsgen += endmember_line + "\n"

    non_endmember_entries = [(level_num, vals) for level_num, vals in comp_entries if level_num > 0 and not _is_pure_one(vals)]

    written: set[str] = set()

    for level_num, vals in non_endmember_entries:
        if level_num > level:
            continue

        for active_letter in all_letters:
            line = f"level={level_num}"

            for letter in all_letters:
                if letter == active_letter:
                    line += f"\t\t{letter}={_fmt_comp(vals)}"
                else:
                    line += f"\t\t{letter}=1.0"

            if line not in written:
                sqsgen += line + "\n"
                written.add(line)

    for k in range(2, len(all_letters) + 1):
        for active_letters in combinations(all_letters, k):
            for combo in product(non_endmember_entries, repeat=k):
                combo_levels = [entry[0] for entry in combo]
                combo_vals = [entry[1] for entry in combo]

                combined_level = max(combo_levels)

                if combined_level > level:
                    continue

                line = f"level={combined_level}"

                for letter in all_letters:
                    if letter in active_letters:
                        idx = active_letters.index(letter)
                        line += f"\t\t{letter}={_fmt_comp(combo_vals[idx])}"
                    else:
                        line += f"\t\t{letter}=1.0"

                if line not in written:
                    sqsgen += line + "\n"
                    written.add(line)

    # Append fixed-composition sublattices to every generated line
    if fixed_compositions:
        fixed_suffix = "".join(f"\t\t{letter}={_fmt_comp(comp)}" for letter, comp in fixed_compositions.items())
        sqsgen = "\n".join(line + fixed_suffix if line.strip() else line for line in sqsgen.splitlines()) + "\n"

    return sqsgen
