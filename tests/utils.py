
def assert_equations(given, expected, diff=False):

    def clean_str(obj):
        return obj.__str__().replace(' ', '').strip().split('\n')

    given_lines = sorted(clean_str(given))
    expected_lines = sorted(clean_str(expected))

    if diff:
        print("diff:")
        for g, e in zip(given_lines, expected_lines):
            print('given:    {}'.format(g))
            print('expected: {}'.format(e))

    assert given_lines == expected_lines
