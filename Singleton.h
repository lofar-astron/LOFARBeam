template<typename T>
class Singleton
{
    public:
        static std::shared_ptr<T> getInstance()
        {
            static std::shared_ptr<T>  instance(new T()); // Guaranteed to be destroyed.
                                                          // Instantiated on first use.
            return instance;
        }
    private:
        Singleton() {}                    // Constructor? (the {} brackets) are needed here.

    public:
        Singleton(Singleton const&)       = delete;
        void operator=(Singleton const&)  = delete;

};
