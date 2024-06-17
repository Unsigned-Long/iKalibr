/*! \file yaml.hpp
    \brief YAML input and output archives */
#ifndef IKALIBR_CEREAL_YAML_HPP_
#define IKALIBR_CEREAL_YAML_HPP_

#include <cereal/cereal.hpp>
#include <cereal/details/util.hpp>
#include <cereal/external/base64.hpp>

#include <yaml-cpp/yaml.h>

#include <limits>
#include <sstream>
#include <stack>
#include <vector>
#include <string>
#include <cstring>

namespace cereal {

class YAMLOutputArchive : public OutputArchive<YAMLOutputArchive>, public traits::TextArchive {
    enum class NodeType { StartObject, InObject, StartArray, InArray };

public:
    YAMLOutputArchive(std::ostream &stream)
        : OutputArchive<YAMLOutputArchive>(this),
          out(stream),
          nextName(nullptr) {
        nameCounter.push(0);
        nodeStack.push(NodeType::StartObject);
    }

    ~YAMLOutputArchive() CEREAL_NOEXCEPT {
        if (nodeStack.top() == NodeType::InObject) {
            emitter << YAML::EndMap;
        } else if (nodeStack.top() == NodeType::InArray) {
            emitter << YAML::EndSeq;
        }

        out << emitter.c_str();
    }

    //! Saves some binary data, encoded as a base64 string, with an optional name
    /*! This will create a new node, optionally named, and insert a value that consists of
        the data encoded as a base64 string */
    void saveBinaryValue(const void *data, size_t size, const char *name = nullptr) {
        setNextName(name);
        writeName();

        auto base64string = base64::encode(reinterpret_cast<const unsigned char *>(data), size);
        saveValue(base64string);
    };

    //! @}
    /*! @name Internal Functionality
        Functionality designed for use by those requiring control over the inner mechanisms of
        the YAMLOutputArchive */
    //! @{

    //! Starts a new node in the YAML output
    /*! The node can optionally be given a name by calling setNextName prior
        to creating the node

        Nodes only need to be started for types that are themselves objects or arrays */
    void startNode() {
        writeName();
        nodeStack.push(NodeType::StartObject);
        nameCounter.push(0);
    }

    //! Designates the most recently added node as finished
    void finishNode() {
        // if we ended up serializing an empty object or array, writeName
        // will never have been called - so start and then immediately end
        // the object/array.
        //
        // We'll also end any object/arrays we happen to be in
        switch (nodeStack.top()) {
            case NodeType::StartArray:
                emitter << YAML::BeginSeq;
            case NodeType::InArray:
                emitter << YAML::EndSeq;
                break;
            case NodeType::StartObject:
                emitter << YAML::BeginMap;
            case NodeType::InObject:
                emitter << YAML::EndMap;
                break;
        }

        nodeStack.pop();
        nameCounter.pop();
    }

    //! Sets the name for the next node created with startNode
    void setNextName(const char *name) { nextName = name; }

    //! Saves a bool to the current node
    void saveValue(bool b) { emitter << b; }
    //! Saves an int to the current node
    void saveValue(int i) { emitter << i; }
    //! Saves a uint to the current node
    void saveValue(unsigned u) { emitter << u; }
    //! Saves an int64 to the current node
    void saveValue(int64_t i64) { emitter << i64; }
    //! Saves a uint64 to the current node
    void saveValue(uint64_t u64) { emitter << u64; }
    //! Saves a double to the current node
    void saveValue(double d) { emitter << d; }
    //! Saves a string to the current node
    void saveValue(std::string const &s) { emitter << s; }
    //! Saves a const char * to the current node
    void saveValue(char const *s) { emitter << s; }
    //! Saves a nullptr to the current node
    void saveValue(std::nullptr_t) { emitter << nullptr; }

private:
    // Some compilers/OS have difficulty disambiguating the above for various flavors of longs, so
    // we provide special overloads to handle these cases.

    //! 32 bit signed long saving to current node
    template <class T,
              traits::EnableIf<sizeof(T) == sizeof(std::int32_t), std::is_signed<T>::value> =
                  traits::sfinae>
    inline void saveLong(T l) {
        saveValue(static_cast<std::int32_t>(l));
    }

    //! non 32 bit signed long saving to current node
    template <class T,
              traits::EnableIf<sizeof(T) != sizeof(std::int32_t), std::is_signed<T>::value> =
                  traits::sfinae>
    inline void saveLong(T l) {
        saveValue(static_cast<std::int64_t>(l));
    }

    //! 32 bit unsigned long saving to current node
    template <class T,
              traits::EnableIf<sizeof(T) == sizeof(std::int32_t), std::is_unsigned<T>::value> =
                  traits::sfinae>
    inline void saveLong(T lu) {
        saveValue(static_cast<std::uint32_t>(lu));
    }

    //! non 32 bit unsigned long saving to current node
    template <class T,
              traits::EnableIf<sizeof(T) != sizeof(std::int32_t), std::is_unsigned<T>::value> =
                  traits::sfinae>
    inline void saveLong(T lu) {
        saveValue(static_cast<std::uint64_t>(lu));
    }

public:
#ifdef _MSC_VER
    //! MSVC only long overload to current node
    void saveValue(unsigned long lu) { saveLong(lu); };
#else   // _MSC_VER
    //! Serialize a long if it would not be caught otherwise
    template <class T,
              traits::EnableIf<std::is_same<T, long>::value,
                               !std::is_same<T, std::int32_t>::value,
                               !std::is_same<T, std::int64_t>::value> = traits::sfinae>
    inline void saveValue(T t) {
        saveLong(t);
    }

    //! Serialize an unsigned long if it would not be caught otherwise
    template <class T,
              traits::EnableIf<std::is_same<T, unsigned long>::value,
                               !std::is_same<T, std::uint32_t>::value,
                               !std::is_same<T, std::uint64_t>::value> = traits::sfinae>
    inline void saveValue(T t) {
        saveLong(t);
    }
#endif  // _MSC_VER

    //! Save exotic arithmetic as strings to current node
    /*! Handles long long (if distinct from other types), unsigned long (if distinct), and long
     * double */
    template <class T,
              traits::EnableIf<std::is_arithmetic<T>::value,
                               !std::is_same<T, long>::value,
                               !std::is_same<T, unsigned long>::value,
                               !std::is_same<T, std::int64_t>::value,
                               !std::is_same<T, std::uint64_t>::value,
                               (sizeof(T) >= sizeof(long double) ||
                                sizeof(T) >= sizeof(long long))> = traits::sfinae>
    inline void saveValue(T const &t) {
        std::stringstream ss;
        ss.precision(std::numeric_limits<long double>::max_digits10);
        ss << t;
        saveValue(ss.str());
    }

    //! Write the name of the upcoming node and prepare object/array state
    /*! Since writeName is called for every value that is output, regardless of
        whether it has a name or not, it is the place where we will do a deferred
        check of our node state and decide whether we are in an array or an object.

        The general workflow of saving to the YAML archive is:

        1. (optional) Set the name for the next node to be created, usually done by an NVP
        2. Start the node
        3. (if there is data to save) Write the name of the node (this function)
        4. (if there is data to save) Save the data (with saveValue)
        5. Finish the node
        */
    void writeName() {
        NodeType const &nodeType = nodeStack.top();

        // Start up either an object or an array, depending on state
        if (nodeType == NodeType::StartArray) {
            emitter << YAML::BeginSeq;
            nodeStack.top() = NodeType::InArray;
        } else if (nodeType == NodeType::StartObject) {
            emitter << YAML::BeginMap;
            nodeStack.top() = NodeType::InObject;
        }

        // Array types do not output names
        if (nodeType == NodeType::InArray) {
            return;
        }

        emitter << YAML::Key;

        if (nextName == nullptr) {
            std::string name = "value" + std::to_string(nameCounter.top()++) + "\0";
            saveValue(name);
        } else {
            saveValue(nextName);
            nextName = nullptr;
        }

        emitter << YAML::Value;
    }

    //! Designates that the current node should be output as an array, not an object
    void makeArray() { nodeStack.top() = NodeType::StartArray; }

    //! @}

private:
    std::ostream &out;
    YAML::Emitter emitter;
    char const *nextName;              //!< The next name
    std::stack<uint32_t> nameCounter;  //!< Counter for creating unique names for unnamed nodes
    std::stack<NodeType> nodeStack;

};  // YAMLOutputArchive

class YAMLInputArchive : public InputArchive<YAMLInputArchive>, public traits::TextArchive {
    typedef YAML::const_iterator YAMLIterator;

public:
    YAMLInputArchive(std::istream &stream)
        : InputArchive<YAMLInputArchive>(this),
          itsNextName(nullptr),
          itsDocument(YAML::Load(stream)) {
        if (itsDocument.IsSequence()) {
            itsIteratorStack.emplace_back(itsDocument.begin(), itsDocument.end(),
                                          Iterator::Type::Value);
        }

        else {
            itsIteratorStack.emplace_back(itsDocument.begin(), itsDocument.end(),
                                          Iterator::Type::Member);
        }
    }

    ~YAMLInputArchive() CEREAL_NOEXCEPT = default;

    //! Loads some binary data, encoded as a base64 string
    /*! This will automatically start and finish a node to load the data, and can be called directly
       by users.

        Note that this follows the same ordering rules specified in the class description in regards
        to loading in/out of order */
    void loadBinaryValue(void *data, size_t size, const char *name = nullptr) {
        itsNextName = name;

        std::string encoded;
        loadValue(encoded);
        auto decoded = base64::decode(encoded);

        if (size != decoded.size()) {
            throw Exception("Decoded binary data size does not match specified size");
        }

        std::memcpy(data, decoded.data(), decoded.size());
        itsNextName = nullptr;
    };

private:
    //! @}
    /*! @name Internal Functionality
        Functionality designed for use by those requiring control over the inner mechanisms of
        the YAMLInputArchive */
    //! @{

    //! An internal iterator that handles both array and object types
    class Iterator {
    public:
        enum Type { Value, Member, Null_ };

        Iterator()
            : itsType(Null_) {}

        Iterator(YAMLIterator begin, YAMLIterator end, Type type)
            : itsItBegin(begin),
              itsItEnd(end),
              itsItCurrent(begin),
              itsType(type) {
            if (std::distance(begin, end) == 0) {
                itsType = Null_;
            }

            if (itsType == Member && itsItCurrent != itsItEnd) {
                currentName = itsItCurrent->first.as<std::string>();
            }
        }

        //! Advance to the next node
        Iterator &operator++() {
            ++itsItCurrent;
            if (itsType == Member && itsItCurrent != itsItEnd) {
                currentName = itsItCurrent->first.as<std::string>();
            }
            return *this;
        }

        //! Get the value of the current node
        YAML::Node value() {
            switch (itsType) {
                case Value:
                    return *itsItCurrent;
                case Member:
                    return itsItCurrent->second;
                default:
                    throw cereal::Exception(
                        "YAMLInputArchive internal error: null or empty iterator to object or "
                        "array!");
            }
        }

        //! Get the name of the current node, or nullptr if it has no name
        const char *name() const {
            if (currentName != "") {
                return currentName.c_str();
            } else {
                return nullptr;
            }
        }

        //! Adjust our position such that we are at the node with the given name
        /*! @throws Exception if no such named node exists */
        inline void search(const char *searchName) {
            auto index = 0;
            const auto len = std::strlen(searchName);
            for (itsItCurrent = itsItBegin; itsItCurrent != itsItEnd; ++itsItCurrent, ++index) {
                currentName = itsItCurrent->first.as<std::string>();
                if ((std::strncmp(searchName, currentName.c_str(), len) == 0) &&
                    (std::strlen(currentName.c_str()) == len)) {
                    return;
                }
            }

            throw Exception("YAML Parsing failed - provided NVP (" + std::string(searchName) +
                            ") not found");
        }

    private:
        YAMLIterator itsItBegin, itsItEnd;  //!< Member or value iterator (object/array)
        YAMLIterator itsItCurrent;          //!< Current iterator
        std::string currentName;            //!< Current name
        Type itsType;  //!< Whether this holds values (array) or members (objects) or nothing
    };

    //! Searches for the expectedName node if it doesn't match the actualName
    /*! This needs to be called before every load or node start occurs.  This function will
        check to see if an NVP has been provided (with setNextName) and if so, see if that name
       matches the actual next name given.  If the names do not match, it will search in the current
       level of the JSON for that name. If the name is not found, an exception will be thrown.

        Resets the NVP name after called.

        @throws Exception if an expectedName is given and not found */
    inline void search() {
        // The name an NVP provided with setNextName()
        if (itsNextName) {
            // The actual name of the current node
            auto const actualName = itsIteratorStack.back().name();

            // Do a search if we don't see a name coming up, or if the names don't match
            if (!actualName || std::strcmp(itsNextName, actualName) != 0) {
                itsIteratorStack.back().search(itsNextName);
            }
        }

        itsNextName = nullptr;
    }

public:
    //! Starts a new node, going into its proper iterator
    /*! This places an iterator for the next node to be parsed onto the iterator stack.  If the next
        node is an array, this will be a value iterator, otherwise it will be a member iterator.

        By default our strategy is to start with the document root node and then recursively iterate
       through all children in the order they show up in the document. We don't need to know NVPs to
       do this; we'll just blindly load in the order things appear in.

        If we were given an NVP, we will search for it if it does not match our the name of the next
       node that would normally be loaded.  This functionality is provided by search(). */
    void startNode() {
        search();

        auto value = itsIteratorStack.back().value();

        if (value.IsSequence()) {
            itsIteratorStack.emplace_back(value.begin(), value.end(), Iterator::Type::Value);
        } else {
            itsIteratorStack.emplace_back(value.begin(), value.end(), Iterator::Type::Member);
        }
    }

    //! Finishes the most recently started node
    void finishNode() {
        itsIteratorStack.pop_back();
        ++itsIteratorStack.back();
    }

    //! Retrieves the current node name
    /*! @return nullptr if no name exists */
    const char *getNodeName() const { return itsIteratorStack.back().name(); }

    //! Sets the name for the next node created with startNode
    void setNextName(const char *name) { itsNextName = name; }

    //! Loads a value from the current node - small signed overload
    template <
        class T,
        traits::EnableIf<std::is_signed<T>::value, sizeof(T) < sizeof(int64_t)> = traits::sfinae>
    inline void loadValue(T &val) {
        search();

        val = static_cast<T>(itsIteratorStack.back().value().as<int>());
        ++itsIteratorStack.back();
    }

    //! Loads a value from the current node - small unsigned overload
    template <class T,
              traits::EnableIf<std::is_unsigned<T>::value,
                               sizeof(T) < sizeof(uint64_t),
                               !std::is_same<bool, T>::value> = traits::sfinae>
    inline void loadValue(T &val) {
        search();

        val = static_cast<T>(itsIteratorStack.back().value().as<unsigned int>());
        ++itsIteratorStack.back();
    }

    //! Loads a value from the current node - bool overload
    void loadValue(bool &val) {
        search();
        val = itsIteratorStack.back().value().as<bool>();
        ++itsIteratorStack.back();
    }
    //! Loads a value from the current node - int64 overload
    void loadValue(int64_t &val) {
        search();
        val = itsIteratorStack.back().value().as<int64_t>();
        ++itsIteratorStack.back();
    }
    //! Loads a value from the current node - uint64 overload
    void loadValue(uint64_t &val) {
        search();
        val = itsIteratorStack.back().value().as<uint64_t>();
        ++itsIteratorStack.back();
    }
    //! Loads a value from the current node - float overload
    void loadValue(float &val) {
        search();
        val = static_cast<float>(itsIteratorStack.back().value().as<float>());
        ++itsIteratorStack.back();
    }
    //! Loads a value from the current node - double overload
    void loadValue(double &val) {
        search();
        val = itsIteratorStack.back().value().as<double>();
        ++itsIteratorStack.back();
    }
    //! Loads a value from the current node - string overload
    void loadValue(std::string &val) {
        search();
        val = itsIteratorStack.back().value().as<std::string>();
        ++itsIteratorStack.back();
    }
    //! Loads a nullptr from the current node
    void loadValue(std::nullptr_t &) {
        search(); /*TodoCEREAL_RAPIDJSON_ASSERT(itsIteratorStack.back().value().IsNull());*/
        ++itsIteratorStack.back();
    }

// Special cases to handle various flavors of long, which tend to conflict with
// the int32_t or int64_t on various compiler/OS combinations.  MSVC doesn't need any of this.
#ifndef _MSC_VER
private:
    //! 32 bit signed long loading from current node
    template <class T>
    inline typename std::enable_if<sizeof(T) == sizeof(std::int32_t) && std::is_signed<T>::value,
                                   void>::type
    loadLong(T &l) {
        loadValue(reinterpret_cast<std::int32_t &>(l));
    }

    //! non 32 bit signed long loading from current node
    template <class T>
    inline typename std::enable_if<sizeof(T) == sizeof(std::int64_t) && std::is_signed<T>::value,
                                   void>::type
    loadLong(T &l) {
        loadValue(reinterpret_cast<std::int64_t &>(l));
    }

    //! 32 bit unsigned long loading from current node
    template <class T>
    inline typename std::enable_if<sizeof(T) == sizeof(std::uint32_t) && !std::is_signed<T>::value,
                                   void>::type
    loadLong(T &lu) {
        loadValue(reinterpret_cast<std::uint32_t &>(lu));
    }

    //! non 32 bit unsigned long loading from current node
    template <class T>
    inline typename std::enable_if<sizeof(T) == sizeof(std::uint64_t) && !std::is_signed<T>::value,
                                   void>::type
    loadLong(T &lu) {
        loadValue(reinterpret_cast<std::uint64_t &>(lu));
    }

public:
    //! Serialize a long if it would not be caught otherwise
    template <class T>
    inline
        typename std::enable_if<std::is_same<T, long>::value && sizeof(T) >= sizeof(std::int64_t) &&
                                    !std::is_same<T, std::int64_t>::value,
                                void>::type
        loadValue(T &t) {
        loadLong(t);
    }

    //! Serialize an unsigned long if it would not be caught otherwise
    template <class T>
    inline typename std::enable_if<std::is_same<T, unsigned long>::value &&
                                       sizeof(T) >= sizeof(std::uint64_t) &&
                                       !std::is_same<T, std::uint64_t>::value,
                                   void>::type
    loadValue(T &t) {
        loadLong(t);
    }
#endif  // _MSC_VER

private:
    //! Convert a string to a long long
    void stringToNumber(std::string const &str, long long &val) { val = std::stoll(str); }
    //! Convert a string to an unsigned long long
    void stringToNumber(std::string const &str, unsigned long long &val) { val = std::stoull(str); }
    //! Convert a string to a long double
    void stringToNumber(std::string const &str, long double &val) { val = std::stold(str); }

public:
    //! Loads a value from the current node - long double and long long overloads
    template <class T,
              traits::EnableIf<std::is_arithmetic<T>::value,
                               !std::is_same<T, long>::value,
                               !std::is_same<T, unsigned long>::value,
                               !std::is_same<T, std::int64_t>::value,
                               !std::is_same<T, std::uint64_t>::value,
                               (sizeof(T) >= sizeof(long double) ||
                                sizeof(T) >= sizeof(long long))> = traits::sfinae>
    inline void loadValue(T &val) {
        std::string encoded;
        loadValue(encoded);
        stringToNumber(encoded, val);
    }

    //! Loads the size for a SizeTag
    void loadSize(size_type &size) {
        if (itsIteratorStack.size() == 1) {
            size = itsDocument.size();
        } else {
            size = (itsIteratorStack.rbegin() + 1)->value().size();
        }
    }

    //! @}

private:
    const char *itsNextName;                 //!< Next name set by NVP
    std::vector<Iterator> itsIteratorStack;  //!< 'Stack' of YAML iterators
    YAML::Node itsDocument;                  //!< YAML document

};  // YAMLInputArchive

// ######################################################################
// YAMLArchive prologue and epilogue functions
// ######################################################################

// ######################################################################
//! Prologue for NVPs for YAML archives
/*! NVPs do not start or finish nodes - they just set up the names */
template <class T>
inline void prologue(YAMLOutputArchive &, NameValuePair<T> const &) {}

//! Prologue for NVPs for YAML archives
template <class T>
inline void prologue(YAMLInputArchive &, NameValuePair<T> const &) {}

// ######################################################################
//! Epilogue for NVPs for YAML archives
/*! NVPs do not start or finish nodes - they just set up the names */
template <class T>
inline void epilogue(YAMLOutputArchive &, NameValuePair<T> const &) {}

//! Epilogue for NVPs for YAML archives
/*! NVPs do not start or finish nodes - they just set up the names */
template <class T>
inline void epilogue(YAMLInputArchive &, NameValuePair<T> const &) {}

// ######################################################################
//! Prologue for SizeTags for YAML archives
/*! SizeTags are strictly ignored for YAML, they just indicate
    that the current node should be made into an array */
template <class T>
inline void prologue(YAMLOutputArchive &ar, SizeTag<T> const &) {
    ar.makeArray();
}

//! Prologue for SizeTags for YAML archives
template <class T>
inline void prologue(YAMLInputArchive &, SizeTag<T> const &) {}

// ######################################################################
//! Epilogue for SizeTags for YAML archives
/*! SizeTags are strictly ignored for YAML */
template <class T>
inline void epilogue(YAMLOutputArchive &, SizeTag<T> const &) {}

//! Epilogue for SizeTags for YAML archives
template <class T>
inline void epilogue(YAMLInputArchive &, SizeTag<T> const &) {}

// ######################################################################
//! Prologue for all other types for YAML archives (except minimal types)
/*! Starts a new node, named either automatically or by some NVP,
    that may be given data by the type about to be archived

    Minimal types do not start or finish nodes */
template <
    class T,
    traits::EnableIf<
        !std::is_arithmetic<T>::value,
        !traits::has_minimal_base_class_serialization<T,
                                                      traits::has_minimal_output_serialization,
                                                      YAMLOutputArchive>::value,
        !traits::has_minimal_output_serialization<T, YAMLOutputArchive>::value> = traits::sfinae>
inline void prologue(YAMLOutputArchive &ar, T const &) {
    ar.startNode();
}

//! Prologue for all other types for YAML archives
template <
    class T,
    traits::EnableIf<
        !std::is_arithmetic<T>::value,
        !traits::has_minimal_base_class_serialization<T,
                                                      traits::has_minimal_input_serialization,
                                                      YAMLInputArchive>::value,
        !traits::has_minimal_input_serialization<T, YAMLInputArchive>::value> = traits::sfinae>
inline void prologue(YAMLInputArchive &ar, T const &) {
    ar.startNode();
}

// ######################################################################
//! Epilogue for all other types other for YAML archives (except minimal types)
/*! Finishes the node created in the prologue

    Minimal types do not start or finish nodes */
template <
    class T,
    traits::EnableIf<
        !std::is_arithmetic<T>::value,
        !traits::has_minimal_base_class_serialization<T,
                                                      traits::has_minimal_output_serialization,
                                                      YAMLOutputArchive>::value,
        !traits::has_minimal_output_serialization<T, YAMLOutputArchive>::value> = traits::sfinae>
inline void epilogue(YAMLOutputArchive &ar, T const &) {
    ar.finishNode();
}

//! Epilogue for all other types other for YAML archives
template <
    class T,
    traits::EnableIf<
        !std::is_arithmetic<T>::value,
        !traits::has_minimal_base_class_serialization<T,
                                                      traits::has_minimal_input_serialization,
                                                      YAMLInputArchive>::value,
        !traits::has_minimal_input_serialization<T, YAMLInputArchive>::value> = traits::sfinae>
inline void epilogue(YAMLInputArchive &ar, T const &) {
    ar.finishNode();
}

// ######################################################################
//! Prologue for arithmetic types for YAML archives
inline void prologue(YAMLOutputArchive &ar, std::nullptr_t const &) { ar.writeName(); }

//! Prologue for arithmetic types for YAML archives
inline void prologue(YAMLInputArchive &, std::nullptr_t const &) {}

// ######################################################################
//! Epilogue for arithmetic types for YAML archives
inline void epilogue(YAMLOutputArchive &, std::nullptr_t const &) {}

//! Epilogue for arithmetic types for YAML archives
inline void epilogue(YAMLInputArchive &, std::nullptr_t const &) {}

// ######################################################################
//! Prologue for arithmetic types for YAML archives
template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae>
inline void prologue(YAMLOutputArchive &ar, T const &) {
    ar.writeName();
}

//! Prologue for arithmetic types for YAML archives
template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae>
inline void prologue(YAMLInputArchive &, T const &) {}

// ######################################################################
//! Epilogue for arithmetic types for YAML archives
template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae>
inline void epilogue(YAMLOutputArchive &, T const &) {}

//! Epilogue for arithmetic types for YAML archives
template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae>
inline void epilogue(YAMLInputArchive &, T const &) {}

// ######################################################################
//! Prologue for strings for YAML archives
template <class CharT, class Traits, class Alloc>
inline void prologue(YAMLOutputArchive &ar, std::basic_string<CharT, Traits, Alloc> const &) {
    ar.writeName();
}

//! Prologue for strings for YAML archives
template <class CharT, class Traits, class Alloc>
inline void prologue(YAMLInputArchive &, std::basic_string<CharT, Traits, Alloc> const &) {}

// ######################################################################
//! Epilogue for strings for YAML archives
template <class CharT, class Traits, class Alloc>
inline void epilogue(YAMLOutputArchive &, std::basic_string<CharT, Traits, Alloc> const &) {}

//! Epilogue for strings for YAML archives
template <class CharT, class Traits, class Alloc>
inline void epilogue(YAMLInputArchive &, std::basic_string<CharT, Traits, Alloc> const &) {}

// ######################################################################
// Common YAMLArchive serialization functions
// ######################################################################
//! Serializing NVP types to YAML
template <class T>
inline void CEREAL_SAVE_FUNCTION_NAME(YAMLOutputArchive &ar, NameValuePair<T> const &t) {
    ar.setNextName(t.name);
    ar(t.value);
}

template <class T>
inline void CEREAL_LOAD_FUNCTION_NAME(YAMLInputArchive &ar, NameValuePair<T> &t) {
    ar.setNextName(t.name);
    ar(t.value);
}

//! Saving for nullptr to YAML
inline void CEREAL_SAVE_FUNCTION_NAME(YAMLOutputArchive &ar, std::nullptr_t const &t) {
    ar.saveValue(t);
}

//! Loading arithmetic from YAML
inline void CEREAL_LOAD_FUNCTION_NAME(YAMLInputArchive &ar, std::nullptr_t &t) { ar.loadValue(t); }

//! Saving for arithmetic to YAML
template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae>
inline void CEREAL_SAVE_FUNCTION_NAME(YAMLOutputArchive &ar, T const &t) {
    ar.saveValue(t);
}

//! Loading arithmetic from YAML
template <class T, traits::EnableIf<std::is_arithmetic<T>::value> = traits::sfinae>
inline void CEREAL_LOAD_FUNCTION_NAME(YAMLInputArchive &ar, T &t) {
    ar.loadValue(t);
}

//! saving string to YAML
template <class CharT, class Traits, class Alloc>
inline void CEREAL_SAVE_FUNCTION_NAME(YAMLOutputArchive &ar,
                                      std::basic_string<CharT, Traits, Alloc> const &str) {
    ar.saveValue(str);
}

//! loading string from YAML
template <class CharT, class Traits, class Alloc>
inline void CEREAL_LOAD_FUNCTION_NAME(YAMLInputArchive &ar,
                                      std::basic_string<CharT, Traits, Alloc> &str) {
    ar.loadValue(str);
}

// ######################################################################
//! Saving SizeTags to YAML
template <class T>
inline void CEREAL_SAVE_FUNCTION_NAME(YAMLOutputArchive &, SizeTag<T> const &) {
    // nothing to do here, we don't explicitly save the size
}

//! Loading SizeTags from YAML
template <class T>
inline void CEREAL_LOAD_FUNCTION_NAME(YAMLInputArchive &ar, SizeTag<T> &st) {
    ar.loadSize(st.size);
}

}  // namespace cereal

// register archives for polymorphic support
CEREAL_REGISTER_ARCHIVE(cereal::YAMLInputArchive)
CEREAL_REGISTER_ARCHIVE(cereal::YAMLOutputArchive)

// tie input and output archives together
CEREAL_SETUP_ARCHIVE_TRAITS(cereal::YAMLInputArchive, cereal::YAMLOutputArchive)

#endif  // IKALIBR_CEREAL_YAML_HPP_
