
def assert_equations(given, expected):

    def clean_str(obj):
        return obj.__str__().replace(' ', '').strip().split('\n')

    given_lines = sorted(clean_str(given))
    expected_lines = sorted(clean_str(expected))

    diffs = []
    for g, e in zip(given_lines, expected_lines):
        if g == e:
            diffs.append('correct:\t{}'.format(g))
        else:
            diffs.append('given:\t\t{}'.format(g))
            diffs.append('expected:\t{}'.format(e))
    diffs = '\n'.join(diffs)

    assert given_lines == expected_lines, 'equations are different\n{}'.format(diffs)
