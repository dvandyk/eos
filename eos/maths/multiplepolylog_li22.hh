// This file is made by Hjalte Frellesvig, Damiano Tommasini, and Christopher Wever
// For further details see ArXiv:1601.02649
// Last modification: 28th Jan. 2016

namespace eos
{
    namespace li22_impl
    {
        // BELOW ARE THE COEFFICIENTS FOR LI22FAST

        // cc1 = Table[(i + 1)^2/((i + 2)^2), {i, 1, 100}];
        // cc2 = Table[1/(i + 1)^2/((i + 2)^2), {i, 1, 100}];

        extern const double ccli221[100];
        extern const double ccli222[100];

        // FIXME

        extern const double ccli2logA1[100];
        extern const double ccli3logA1[100];
        extern const double ccli4logA1[100];
    }
}