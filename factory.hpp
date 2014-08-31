#ifndef FF_FLAME_FACTORY_HPP
#define FF_FLAME_FACTORY_HPP

#include <map>

template < class VariantType,
	       typename KeyType,
	       typename VariantCreator> 
class flame_factory
{
public:
	VariantType* create_variant(const KeyType& id)
    {
        auto variant_iter = creator_map.find(id);
        if(variant_iter != creator_map.end())
            return (variant_iter->second)();
        return nullptr;
    }

	bool register_variant(const KeyType& id, VariantCreator creator)
    {
        return creator_map.insert(std::pair<KeyType, VariantCreator>(id, creator)).second;
    }

	bool unregister_variant(const KeyType& id)
    {
        return (creator_map.erase(id) == 1);
    }

private:
	std::map<KeyType, VariantCreator> creator_map;
};

#endif
