// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_
#define SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_

#include <seqan/map.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <seqan/arg_parse/arg_parse_type_support.h>
#include <seqan/arg_parse/arg_parse_argument.h>
#include <seqan/arg_parse/arg_parse_option.h>

#include <seqan/misc/misc_terminal.h>
#include <seqan/misc/tool_doc.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// friend declaration to make addOption() and hideOption() available
// in ArgumentParser::init()
class ArgumentParser;
class ArgParseOption;
void addOption(ArgumentParser & me, ArgParseOption const & opt);
void hideOption(ArgumentParser & me, std::string const & name, bool hide);

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

/**
.Class.ArgumentParser
..cat:Miscellaneous
..summary:Stores multiple @Class.ArgParseOption@ objects and parses the command line arguments for these options.
..signature:ArgumentParser
..include:seqan/arg_parse.h
..remarks:
See the documentation of @Class.ToolDoc@ on how to format text.
Where possible, formatting is added automatically for you.
You have to use formatting in the following places: (1) usage lines, (2) option help texts, (3) description and additional text sections.
..example.text:
The following gives a simple example of how to use the @Class.ArgumentParser@.
..example.code:
ArgumentParser parser("alf");
setShortDescription(parser, "Alignment free sequence comparison");
setVersion(parser, "1.0");
setDate(parser, "Jan 2010");

addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN\\fP \\fB-o\\fP \\fIOUT\\fP");

addDescription(parser,
               "ALF can be used to calculate the pairwise similarity of sequences "
               "using alignment-free methods. All methods which are implemented are "
               "based on k-mer counts.");

addOption(parser, ArgParseOption("i", "inputFile", "Name of the multi-FASTA input.",
                                 ArgParseArgument(ArgParseArgument::INPUTFILE, false, "IN")));
setRequired(parser, "i");

addOption(parser, ArgParseOption("o", "outputFile", "Name of the multi-FASTA input.",
                                 ArgParseArgument(ArgParseArgument::OUTPUTFILE, false, "OUT")));
setRequired(parser, "o");

addTextSection(parser, "See Also");
addText(parser, "http://www.seqan.de/projects/alf");
..see:Class.ToolDoc

.Memfunc.ArgumentParser#ArgumentParser
..class:Class.ArgumentParser
..summary:Constructor
..signature:ArgumentParser ()
..signature:ArgumentParser (applicationName)
..param.applicationName:A std::string containing the name of the application.
..remarks:If the name of the application is not passed to the constructor it will be extracted from the command line.
*/

class ArgumentParser
{
public:

    // ----------------------------------------------------------------------------
    // Enum ParseResult
    // ----------------------------------------------------------------------------

    // will be used as return value of parse(..) to indicate whether parsing worked
    enum ParseResult
    {
        PARSE_OK,
        PARSE_ERROR,
        PARSE_HELP,
        PARSE_VERSION,
        PARSE_WRITE_CTD,
        PARSE_EXPORT_HELP
    };

    // ----------------------------------------------------------------------------
    // Class Typedefs
    // ----------------------------------------------------------------------------

    typedef std::vector<ArgParseOption>   TOptionMap;
    typedef std::vector<ArgParseArgument> TArgumentMap;
    typedef Size<TOptionMap>::Type        TOptionMapSize;
    typedef Size<TArgumentMap>::Type      TArgumentMapSize;

    typedef std::map<std::string, TOptionMapSize> TStringMap;
    typedef std::vector<std::string>              TValueMap;

    // ----------------------------------------------------------------------------
    // Mapping of option names to options
    // ----------------------------------------------------------------------------

    TStringMap   shortNameMap;
    TStringMap   longNameMap;
    TOptionMap   optionMap;
    TArgumentMap argumentList;

    // ----------------------------------------------------------------------------
    // Documentation Members
    // ----------------------------------------------------------------------------

    ToolDoc                  _toolDoc;      // the tool doc for all user specified
                                            // text
    ToolDoc                  _description;  // the description which we need to
                                            // separate to put it on top of the rest
    std::vector<std::string> _usageText;    // the usage lines as strings, to avoid
                                            // interference with the rest of the doc

    // ----------------------------------------------------------------------------
    // Function init()
    // ----------------------------------------------------------------------------

    void init()
    {
        addOption(*this, ArgParseOption("h", "help", "Displays this help message."));

        // hidden flags used for export of man pages and ctd formats
        addOption(*this, ArgParseOption("", "write-ctd", "Exports the app's interface description to a .ctd file.", ArgParseArgument::OUTPUTFILE));
        hideOption(*this, "write-ctd", true);
        addOption(*this, ArgParseOption("", "export-help", "Export help to a format. One of {'html', 'man', 'txt'}.", ArgParseArgument::STRING, false, "FORMAT"));
        hideOption(*this, "export-help", true);

        // this is our ToolDoc only for the Description, we will later append it to the
        // real ToolDoc, but we need to separate it to ease the formating
        addSection(_description, "Description");
    }

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    ArgumentParser()
    {
        init();
    }

    ArgumentParser(std::string const & _appName)
    {
        setName(_toolDoc, _appName);
        init();
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function hasOption()
// ----------------------------------------------------------------------------

/**
.Function.hasOption:
..summary:Returns whether a certain option is registered in the parser.
..cat:Miscellaneous
..signature:hasOption(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the option.
..returns:$true$ if the option is registered.
..include:seqan/arg_parse.h
*/

inline bool hasOption(ArgumentParser const & me, std::string const & name)
{
    return hasKey(me.shortNameMap, name) || hasKey(me.longNameMap, name);
}

// ----------------------------------------------------------------------------
// Function addOption()
// ----------------------------------------------------------------------------

/**
.Function.addOption
..summary:Adds a @Class.ArgParseOption@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addOption(parser, option)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The new @Class.ArgParseOption@ object that should be added.
...type:Class.ArgParseOption
..include:seqan/arg_parse.h
*/

inline void addOption(ArgumentParser & me, ArgParseOption const & opt)
{
    // check if an option with the same identifiers was already registered
    SEQAN_CHECK(!hasOption(me, opt.shortName), "There already is an option with the name %s!", toCString(opt.shortName));
    SEQAN_CHECK(!hasOption(me, opt.longName), "There already is an option with the name %s!", toCString(opt.longName));

    // finally append the option
    appendValue(me.optionMap, opt);

    if (!empty(opt.shortName))
        me.shortNameMap.insert(std::make_pair<std::string, ArgumentParser::TOptionMapSize>(opt.shortName, length(me.optionMap) - 1));
    if (!empty(opt.longName))
        me.longNameMap.insert(std::make_pair<std::string, ArgumentParser::TOptionMapSize>(opt.longName, length(me.optionMap) - 1));
}

// ----------------------------------------------------------------------------
// Function addArgument()
// ----------------------------------------------------------------------------

/**
.Function.addArgument
..summary:Adds a @Class.ArgParseArgument@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addArgument(parser, argument)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The new @Class.ArgParseArgument@ object that should be added.
...type:Class.ArgParseArgument
..include:seqan/arg_parse.h
*/

inline void addArgument(ArgumentParser & me, ArgParseArgument const & arg)
{
    // check previous arguments
    //  .. lists can only be last argument
    if (!me.argumentList.empty())
    {
        SEQAN_CHECK(!isListArgument(me.argumentList[me.argumentList.size() - 1]),
                    "You cannot add an additional argument after a list argument.");
    }

    // check current argument
    //  .. arguments should not have default values
    SEQAN_CHECK(arg.defaultValue.empty(), "Arguments cannot have default values.");
    SEQAN_CHECK(arg._numberOfValues == 1, "n-Tuple of arguments are not supported.");

    me.argumentList.push_back(arg);
}

// ----------------------------------------------------------------------------
// Function _getOptionIndex()
// ----------------------------------------------------------------------------
// note that it is assumed that the option exists if this method is called

inline ArgumentParser::TOptionMapSize _getOptionIndex(ArgumentParser const & me,
                                                      std::string const & name)
{
    ArgumentParser::TOptionMapSize option_index;
    if (me.shortNameMap.find(name) != me.shortNameMap.end())
    {
        option_index = me.shortNameMap.find(name)->second;
    }
    else
    {
        option_index = me.longNameMap.find(name)->second;
    }
    return option_index;
}

// ----------------------------------------------------------------------------
// Function getOption()
// ----------------------------------------------------------------------------

/**
.Function.getOption
..summary:Returns a reference to the specified option.
..cat:Miscellaneous
..signature:getOption(parser, optionName)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..returns: a reference to the specified @Class.ArgParseOption@ object.
..include:seqan/arg_parse.h
*/

inline ArgParseOption & getOption(ArgumentParser & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return me.optionMap[_getOptionIndex(me, name)];
}

inline ArgParseOption const & getOption(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return me.optionMap[_getOptionIndex(me, name)];
}

// ----------------------------------------------------------------------------
// Function setRequired()
// ----------------------------------------------------------------------------

/**
.Function.setRequired
..summary:Sets whether or not the option defined by the parameter $name$ (which can be
 either the short or the long name) is mandatory.
..cat:Miscellaneous
..signature:setRequired(parser, optionName [, required])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.required:The new required value of the option. Default is true.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void setRequired(ArgumentParser & me, std::string const & name, bool required = true)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return setRequired(getOption(me, name), required);
}

// ----------------------------------------------------------------------------
// Function hideOption()
// ----------------------------------------------------------------------------

/**
.Function.hideOption
..summary:Hides the ArgParseOption defined by the parameter $name$ (which can be
 either the short or the long name) from the help screen.
..cat:Miscellaneous
..signature:hideOption(parser, optionName [, hide])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.hide:The new visibility of the option. Default is false.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void hideOption(ArgumentParser & me, std::string const & name, bool hide = true)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    hideOption(getOption(me, name), hide);
}

// ----------------------------------------------------------------------------
// Function getArgument()
// ----------------------------------------------------------------------------

/**
.Function.getArgument
..summary:Returns a reference to the specified argument.
..cat:Miscellaneous
..signature:getArgument(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: a reference to the specified @Class.ArgParseArgument@ object.
..include:seqan/arg_parse.h
*/

inline ArgParseArgument & getArgument(ArgumentParser & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}

inline ArgParseArgument const & getArgument(ArgumentParser const & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------

/**
.Function.isSet
..summary:Returns whether an option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSet(parser,optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that identifies the option (either short or long name).
..returns:$true$ if the option was set.
..include:seqan/arg_parse.h
*/

inline bool isSet(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return isSet(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function _allRequiredSet()
// ----------------------------------------------------------------------------

inline bool _allRequiredSet(ArgumentParser const & me)
{
    for (unsigned o = 0; o < length(me.optionMap); ++o)
        if (!isSet(me.optionMap[o]) && isRequired(me.optionMap[o]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function _allArgumentsSet()
// -------------------------------------------------------------------------}---

inline bool _allArgumentsSet(ArgumentParser const & me)
{
    for (unsigned a = 0; a < me.argumentList.size(); ++a)
        if (!isSet(me.argumentList[a]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function getOptionValue()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValue:
..summary:Retrieves the value of an option given either the short or long name.
..cat:Miscellaneous
..signature:getOptionValue(value, parser, optionIdentifier[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..param.argNo:If the option is list, the $argNo$-th list element is returned.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool getOptionValue(TValue & val,
                           ArgumentParser const & me,
                           std::string const & name,
                           unsigned argNo)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return _convertArgumentValue(val, getOption(me, name), getArgumentValue(getOption(me, name), argNo));
}

template <typename TValue>
inline bool getOptionValue(TValue & val,
                           ArgumentParser const & me,
                           std::string const & name)
{
    return getOptionValue(val, me, name, 0);
}

// ----------------------------------------------------------------------------
// Function getOptionValueCount()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValueCount:
..summary:Returns the number of values stored in the specified option.
..cat:Miscellaneous
..signature:getOptionValueCount(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..returns: The number of values stored for this option.
..include:seqan/arg_parse.h
*/

inline unsigned getOptionValueCount(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return getArgumentValues(getOption(me, name)).size();
}

// ----------------------------------------------------------------------------
// Function getArgumentValueCount()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValueCount:
..summary:Retunrs the number of values stored in the specified option.
..cat:Miscellaneous
..signature:getArgumentValueCount(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: The number of values stored for the specified argument.
..include:seqan/arg_parse.h
*/// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================


inline unsigned getArgumentValueCount(ArgumentParser const & me, unsigned argumentPosition)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return getArgumentValues(getArgument(me, argumentPosition)).size();
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValue:
..summary:Retrieves the value of an argument given by its position.
..cat:Miscellaneous
..signature:getArgumentValue(value, parser, argumentPosition[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..param.argNo:If the argument is a list, the $argNo$-th list element is returned.
..returns: $true$ if the requested argument is set and has the requested type, $false$ otherwise.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool getArgumentValue(TValue & value,
                             ArgumentParser & me,
                             unsigned argumentPosition,
                             unsigned argNo)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return _convertArgumentValue(value, getArgument(me, argumentPosition), getArgumentValue(getArgument(me, argumentPosition), argNo));
}

template <typename TValue>
inline bool getArgumentValue(TValue & value,
                             ArgumentParser & me,
                             unsigned argumentPosition)
{
    return getArgumentValue(value, me, argumentPosition, 0);
}

// ----------------------------------------------------------------------------
// Function getOptionValues()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValues
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getOptionValues(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..returns: A $String<std::string>$ of option values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const & getOptionValues(ArgumentParser & me,
                                                        std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return getArgumentValues(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValues
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getArgumentValues(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: A $String<std::string>$ of argument values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const & getArgumentValues(ArgumentParser & me,
                                                          unsigned argumentPosition)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return getArgumentValues(getArgument(me, argumentPosition));
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/**
.Function.setMinValue
..signature:setMinValue(parser,optionName,minValue)
..signature:setMinValue(parser,argumentPosition,minValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.minValue:A std::string containing a string representation of the minimum value of the @Class.ArgParseOption@.
*/

inline void setMinValue(ArgumentParser & me,
                        std::string const & name,
                        std::string const & _minValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMinValue(getOption(me, name), _minValue);
}

inline void setMinValue(ArgumentParser & me,
                        unsigned argumentPosition,
                        std::string const & _minValue)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setMinValue(getArgument(me, argumentPosition), _minValue);
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/**
.Function.setMaxValue
..signature:setMaxValue(parser,optionName,maxValue)
..signature:setMaxValue(parser,argumentPosition,minValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.maxValue:A std::string containing a string representation of the maximum value of the @Class.ArgParseOption@.
*/

inline void setMaxValue(ArgumentParser & me,
                        std::string const & name,
                        std::string const & _maxValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMaxValue(getOption(me, name), _maxValue);
}

inline void setMaxValue(ArgumentParser & me,
                        unsigned argumentPosition,
                        std::string const & _minValue)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setMaxValue(getArgument(me, argumentPosition), _minValue);
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/**
.Function.setValidValues
..signature:setValidValues(parser,optionName,values)
..signature:setValidValues(parser,argumentPosition,values)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.values:A $std::string$ containing all valid entries for the option.
Alternatively you can pass a string containing all values separated by spaces.
*/

inline void setValidValues(ArgumentParser & me,
                           std::string const & name,
                           std::vector<std::string> const & values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), values);
}

inline void setValidValues(ArgumentParser & me,
                           std::string const & name,
                           std::string const & values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), values);
}

inline void setValidValues(ArgumentParser & me,
                           unsigned argumentPosition,
                           std::vector<std::string> const & values)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setValidValues(getArgument(me, argumentPosition), values);
}

inline void setValidValues(ArgumentParser & me,
                           unsigned argumentPosition,
                           std::string const & values)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setValidValues(getArgument(me, argumentPosition), values);
}

}  // namespace seqan

#endif // SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_
// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_
#define SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_

#include <seqan/map.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <seqan/arg_parse/arg_parse_type_support.h>
#include <seqan/arg_parse/arg_parse_argument.h>
#include <seqan/arg_parse/arg_parse_option.h>

#include <seqan/misc/misc_terminal.h>
#include <seqan/misc/tool_doc.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// friend declaration to make addOption() and hideOption() available
// in ArgumentParser::init()
class ArgumentParser;
class ArgParseOption;
void addOption(ArgumentParser & me, ArgParseOption const & opt);
void hideOption(ArgumentParser & me, std::string const & name, bool hide);

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

/**
.Class.ArgumentParser
..cat:Miscellaneous
..summary:Stores multiple @Class.ArgParseOption@ objects and parses the command line arguments for these options.
..signature:ArgumentParser
..include:seqan/arg_parse.h
..remarks:
See the documentation of @Class.ToolDoc@ on how to format text.
Where possible, formatting is added automatically for you.
You have to use formatting in the following places: (1) usage lines, (2) option help texts, (3) description and additional text sections.
..example.text:
The following gives a simple example of how to use the @Class.ArgumentParser@.
..example.code:
ArgumentParser parser("alf");
setShortDescription(parser, "Alignment free sequence comparison");
setVersion(parser, "1.0");
setDate(parser, "Jan 2010");

addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN\\fP \\fB-o\\fP \\fIOUT\\fP");

addDescription(parser,
               "ALF can be used to calculate the pairwise similarity of sequences "
               "using alignment-free methods. All methods which are implemented are "
               "based on k-mer counts.");

addOption(parser, ArgParseOption("i", "inputFile", "Name of the multi-FASTA input.",
                                 ArgParseArgument(ArgParseArgument::INPUTFILE, false, "IN")));
setRequired(parser, "i");

addOption(parser, ArgParseOption("o", "outputFile", "Name of the multi-FASTA input.",
                                 ArgParseArgument(ArgParseArgument::OUTPUTFILE, false, "OUT")));
setRequired(parser, "o");

addTextSection(parser, "See Also");
addText(parser, "http://www.seqan.de/projects/alf");
..see:Class.ToolDoc

.Memfunc.ArgumentParser#ArgumentParser
..class:Class.ArgumentParser
..summary:Constructor
..signature:ArgumentParser ()
..signature:ArgumentParser (applicationName)
..param.applicationName:A std::string containing the name of the application.
..remarks:If the name of the application is not passed to the constructor it will be extracted from the command line.
*/

class ArgumentParser
{
public:

    // ----------------------------------------------------------------------------
    // Enum ParseResult
    // ----------------------------------------------------------------------------

    // will be used as return value of parse(..) to indicate whether parsing worked
    enum ParseResult
    {
        PARSE_OK,
        PARSE_ERROR,
        PARSE_HELP,
        PARSE_VERSION,
        PARSE_WRITE_CTD,
        PARSE_EXPORT_HELP
    };

    // ----------------------------------------------------------------------------
    // Class Typedefs
    // ----------------------------------------------------------------------------

    typedef std::vector<ArgParseOption>   TOptionMap;
    typedef std::vector<ArgParseArgument> TArgumentMap;
    typedef Size<TOptionMap>::Type        TOptionMapSize;
    typedef Size<TArgumentMap>::Type      TArgumentMapSize;

    typedef std::map<std::string, TOptionMapSize> TStringMap;
    typedef std::vector<std::string>              TValueMap;

    // ----------------------------------------------------------------------------
    // Mapping of option names to options
    // ----------------------------------------------------------------------------

    TStringMap   shortNameMap;
    TStringMap   longNameMap;
    TOptionMap   optionMap;
    TArgumentMap argumentList;

    // ----------------------------------------------------------------------------
    // Documentation Members
    // ----------------------------------------------------------------------------

    ToolDoc                  _toolDoc;      // the tool doc for all user specified
                                            // text
    ToolDoc                  _description;  // the description which we need to
                                            // separate to put it on top of the rest
    std::vector<std::string> _usageText;    // the usage lines as strings, to avoid
                                            // interference with the rest of the doc

    // ----------------------------------------------------------------------------
    // Function init()
    // ----------------------------------------------------------------------------

    void init()
    {
        addOption(*this, ArgParseOption("h", "help", "Displays this help message."));

        // hidden flags used for export of man pages and ctd formats
        addOption(*this, ArgParseOption("", "write-ctd", "Exports the app's interface description to a .ctd file.", ArgParseArgument::OUTPUTFILE));
        hideOption(*this, "write-ctd", true);
        addOption(*this, ArgParseOption("", "export-help", "Export help to a format. One of {'html', 'man', 'txt'}.", ArgParseArgument::STRING, false, "FORMAT"));
        hideOption(*this, "export-help", true);

        // this is our ToolDoc only for the Description, we will later append it to the
        // real ToolDoc, but we need to separate it to ease the formating
        addSection(_description, "Description");
    }

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    ArgumentParser()
    {
        init();
    }

    ArgumentParser(std::string const & _appName)
    {
        setName(_toolDoc, _appName);
        init();
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function hasOption()
// ----------------------------------------------------------------------------

/**
.Function.hasOption:
..summary:Returns whether a certain option is registered in the parser.
..cat:Miscellaneous
..signature:hasOption(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the option.
..returns:$true$ if the option is registered.
..include:seqan/arg_parse.h
*/

inline bool hasOption(ArgumentParser const & me, std::string const & name)
{
    return hasKey(me.shortNameMap, name) || hasKey(me.longNameMap, name);
}

// ----------------------------------------------------------------------------
// Function addOption()
// ----------------------------------------------------------------------------

/**
.Function.addOption
..summary:Adds a @Class.ArgParseOption@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addOption(parser, option)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The new @Class.ArgParseOption@ object that should be added.
...type:Class.ArgParseOption
..include:seqan/arg_parse.h
*/

inline void addOption(ArgumentParser & me, ArgParseOption const & opt)
{
    // check if an option with the same identifiers was already registered
    SEQAN_CHECK(!hasOption(me, opt.shortName), "There already is an option with the name %s!", toCString(opt.shortName));
    SEQAN_CHECK(!hasOption(me, opt.longName), "There already is an option with the name %s!", toCString(opt.longName));

    // finally append the option
    appendValue(me.optionMap, opt);

    if (!empty(opt.shortName))
        me.shortNameMap.insert(std::make_pair<std::string, ArgumentParser::TOptionMapSize>(opt.shortName, length(me.optionMap) - 1));
    if (!empty(opt.longName))
        me.longNameMap.insert(std::make_pair<std::string, ArgumentParser::TOptionMapSize>(opt.longName, length(me.optionMap) - 1));
}

// ----------------------------------------------------------------------------
// Function addArgument()
// ----------------------------------------------------------------------------

/**
.Function.addArgument
..summary:Adds a @Class.ArgParseArgument@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addArgument(parser, argument)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The new @Class.ArgParseArgument@ object that should be added.
...type:Class.ArgParseArgument
..include:seqan/arg_parse.h
*/

inline void addArgument(ArgumentParser & me, ArgParseArgument const & arg)
{
    // check previous arguments
    //  .. lists can only be last argument
    if (!me.argumentList.empty())
    {
        SEQAN_CHECK(!isListArgument(me.argumentList[me.argumentList.size() - 1]),
                    "You cannot add an additional argument after a list argument.");
    }

    // check current argument
    //  .. arguments should not have default values
    SEQAN_CHECK(arg.defaultValue.empty(), "Arguments cannot have default values.");
    SEQAN_CHECK(arg._numberOfValues == 1, "n-Tuple of arguments are not supported.");

    me.argumentList.push_back(arg);
}

// ----------------------------------------------------------------------------
// Function _getOptionIndex()
// ----------------------------------------------------------------------------
// note that it is assumed that the option exists if this method is called

inline ArgumentParser::TOptionMapSize _getOptionIndex(ArgumentParser const & me,
                                                      std::string const & name)
{
    ArgumentParser::TOptionMapSize option_index;
    if (me.shortNameMap.find(name) != me.shortNameMap.end())
    {
        option_index = me.shortNameMap.find(name)->second;
    }
    else
    {
        option_index = me.longNameMap.find(name)->second;
    }
    return option_index;
}

// ----------------------------------------------------------------------------
// Function getOption()
// ----------------------------------------------------------------------------

/**
.Function.getOption
..summary:Returns a reference to the specified option.
..cat:Miscellaneous
..signature:getOption(parser, optionName)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..returns: a reference to the specified @Class.ArgParseOption@ object.
..include:seqan/arg_parse.h
*/

inline ArgParseOption & getOption(ArgumentParser & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return me.optionMap[_getOptionIndex(me, name)];
}

inline ArgParseOption const & getOption(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return me.optionMap[_getOptionIndex(me, name)];
}

// ----------------------------------------------------------------------------
// Function setRequired()
// ----------------------------------------------------------------------------

/**
.Function.setRequired
..summary:Sets whether or not the option defined by the parameter $name$ (which can be
 either the short or the long name) is mandatory.
..cat:Miscellaneous
..signature:setRequired(parser, optionName [, required])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.required:The new required value of the option. Default is true.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void setRequired(ArgumentParser & me, std::string const & name, bool required = true)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return setRequired(getOption(me, name), required);
}

// ----------------------------------------------------------------------------
// Function hideOption()
// ----------------------------------------------------------------------------

/**
.Function.hideOption
..summary:Hides the ArgParseOption defined by the parameter $name$ (which can be
 either the short or the long name) from the help screen.
..cat:Miscellaneous
..signature:hideOption(parser, optionName [, hide])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.hide:The new visibility of the option. Default is false.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void hideOption(ArgumentParser & me, std::string const & name, bool hide = true)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    hideOption(getOption(me, name), hide);
}

// ----------------------------------------------------------------------------
// Function getArgument()
// ----------------------------------------------------------------------------

/**
.Function.getArgument
..summary:Returns a reference to the specified argument.
..cat:Miscellaneous
..signature:getArgument(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: a reference to the specified @Class.ArgParseArgument@ object.
..include:seqan/arg_parse.h
*/

inline ArgParseArgument & getArgument(ArgumentParser & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}

inline ArgParseArgument const & getArgument(ArgumentParser const & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------

/**
.Function.isSet
..summary:Returns whether an option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSet(parser,optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that identifies the option (either short or long name).
..returns:$true$ if the option was set.
..include:seqan/arg_parse.h
*/

inline bool isSet(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return isSet(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function _allRequiredSet()
// ----------------------------------------------------------------------------

inline bool _allRequiredSet(ArgumentParser const & me)
{
    for (unsigned o = 0; o < length(me.optionMap); ++o)
        if (!isSet(me.optionMap[o]) && isRequired(me.optionMap[o]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function _allArgumentsSet()
// -------------------------------------------------------------------------}---

inline bool _allArgumentsSet(ArgumentParser const & me)
{
    for (unsigned a = 0; a < me.argumentList.size(); ++a)
        if (!isSet(me.argumentList[a]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function getOptionValue()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValue:
..summary:Retrieves the value of an option given either the short or long name.
..cat:Miscellaneous
..signature:getOptionValue(value, parser, optionIdentifier[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..param.argNo:If the option is list, the $argNo$-th list element is returned.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool getOptionValue(TValue & val,
                           ArgumentParser const & me,
                           std::string const & name,
                           unsigned argNo)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return _convertArgumentValue(val, getOption(me, name), getArgumentValue(getOption(me, name), argNo));
}

template <typename TValue>
inline bool getOptionValue(TValue & val,
                           ArgumentParser const & me,
                           std::string const & name)
{
    return getOptionValue(val, me, name, 0);
}

// ----------------------------------------------------------------------------
// Function getOptionValueCount()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValueCount:
..summary:Returns the number of values stored in the specified option.
..cat:Miscellaneous
..signature:getOptionValueCount(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..returns: The number of values stored for this option.
..include:seqan/arg_parse.h
*/

inline unsigned getOptionValueCount(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return getArgumentValues(getOption(me, name)).size();
}

// ----------------------------------------------------------------------------
// Function getArgumentValueCount()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValueCount:
..summary:Retunrs the number of values stored in the specified option.
..cat:Miscellaneous
..signature:getArgumentValueCount(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: The number of values stored for the specified argument.
..include:seqan/arg_parse.h
*/// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================


inline unsigned getArgumentValueCount(ArgumentParser const & me, unsigned argumentPosition)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return getArgumentValues(getArgument(me, argumentPosition)).size();
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValue:
..summary:Retrieves the value of an argument given by its position.
..cat:Miscellaneous
..signature:getArgumentValue(value, parser, argumentPosition[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..param.argNo:If the argument is a list, the $argNo$-th list element is returned.
..returns: $true$ if the requested argument is set and has the requested type, $false$ otherwise.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool getArgumentValue(TValue & value,
                             ArgumentParser & me,
                             unsigned argumentPosition,
                             unsigned argNo)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return _convertArgumentValue(value, getArgument(me, argumentPosition), getArgumentValue(getArgument(me, argumentPosition), argNo));
}

template <typename TValue>
inline bool getArgumentValue(TValue & value,
                             ArgumentParser & me,
                             unsigned argumentPosition)
{
    return getArgumentValue(value, me, argumentPosition, 0);
}

// ----------------------------------------------------------------------------
// Function getOptionValues()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValues
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getOptionValues(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..returns: A $String<std::string>$ of option values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const & getOptionValues(ArgumentParser & me,
                                                        std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return getArgumentValues(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValues
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getArgumentValues(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: A $String<std::string>$ of argument values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const & getArgumentValues(ArgumentParser & me,
                                                          unsigned argumentPosition)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return getArgumentValues(getArgument(me, argumentPosition));
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/**
.Function.setMinValue
..signature:setMinValue(parser,optionName,minValue)
..signature:setMinValue(parser,argumentPosition,minValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.minValue:A std::string containing a string representation of the minimum value of the @Class.ArgParseOption@.
*/

inline void setMinValue(ArgumentParser & me,
                        std::string const & name,
                        std::string const & _minValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMinValue(getOption(me, name), _minValue);
}

inline void setMinValue(ArgumentParser & me,
                        unsigned argumentPosition,
                        std::string const & _minValue)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setMinValue(getArgument(me, argumentPosition), _minValue);
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/**
.Function.setMaxValue
..signature:setMaxValue(parser,optionName,maxValue)
..signature:setMaxValue(parser,argumentPosition,minValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.maxValue:A std::string containing a string representation of the maximum value of the @Class.ArgParseOption@.
*/

inline void setMaxValue(ArgumentParser & me,
                        std::string const & name,
                        std::string const & _maxValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMaxValue(getOption(me, name), _maxValue);
}

inline void setMaxValue(ArgumentParser & me,
                        unsigned argumentPosition,
                        std::string const & _minValue)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setMaxValue(getArgument(me, argumentPosition), _minValue);
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/**
.Function.setValidValues
..signature:setValidValues(parser,optionName,values)
..signature:setValidValues(parser,argumentPosition,values)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.values:A $std::string$ containing all valid entries for the option.
Alternatively you can pass a string containing all values separated by spaces.
*/

inline void setValidValues(ArgumentParser & me,
                           std::string const & name,
                           std::vector<std::string> const & values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), values);
}

inline void setValidValues(ArgumentParser & me,
                           std::string const & name,
                           std::string const & values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), values);
}

inline void setValidValues(ArgumentParser & me,
                           unsigned argumentPosition,
                           std::vector<std::string> const & values)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setValidValues(getArgument(me, argumentPosition), values);
}

inline void setValidValues(ArgumentParser & me,
                           unsigned argumentPosition,
                           std::string const & values)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    setValidValues(getArgument(me, argumentPosition), values);
}

}  // namespace seqan

#endif // SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_