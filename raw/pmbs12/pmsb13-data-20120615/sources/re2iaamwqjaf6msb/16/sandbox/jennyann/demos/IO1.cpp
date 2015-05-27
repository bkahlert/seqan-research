std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
do
{
    char c = '\0';
    int res = streamReadChar(c, file);
    if (streamEof(file))
        break;  // It's over! But no error.
    if (res)
        return res;  // Pass error code to caller.
}
while (true);