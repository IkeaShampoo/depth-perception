#pragma once
#include <cmath>
#include <iostream>

namespace util {
	template <typename T>
	class linked_list {
		struct node {
			T value;
			node* next;
			node() = default;
			node(T newValue) {
				value = newValue;
				next = 0;
			}
			node* operator [](int i) {
				node* currentNode = this;
				while (i--) { currentNode = currentNode->next; }
				return currentNode;
			}
		};

		void delete_node(node* target) {
			node* nextTarget;
			while (target) {
				nextTarget = target->next;
				delete target;
				target = nextTarget;
			}
		}

		node* node0;
		node* nodeF;
		int length;

		linked_list(node* newNode0, node* newNodeF, int newLength) {
			node0 = newNode0;
			nodeF = newNodeF;
			length = newLength;
		}

		public:

			typedef node* position;

			linked_list(T newValue) {
				node0 = new node(newValue);
				nodeF = node0;
				length = 1;
			}

			linked_list(const T* source, int sourceLength) {
				node0 = new node(source[0]);
				nodeF = node0;
				for (int i = 1; i < sourceLength; i++) {
					nodeF = (nodeF->next = new node(source[i]));
				}
				length = sourceLength;
			}

			linked_list(const linked_list& source) {
				node0 = source.node0;
				nodeF = source.nodeF;
				length = source.length;
			}

			void* operator ()(int i) {
				return (*node0)[i];
			}

			void* getStart() {
				return node0;
			}

			void* next(void* pos) {
				return reinterpret_cast<node*>(pos)->next;
			}

			linked_list operator ()(int start, int end) { // end not included
				node* currentNode = (*node0)[start];
				node* newNode0 = new node(currentNode->value);
				node* newNodeF = newNode0;
				for (int i = start + 1; i < end; i++) {
					currentNode = currentNode->next;
					newNodeF = (newNode0->next = new node(currentNode->value));
				}
				return linked_list(newNode0, newNodeF, end - start);
			}

			linked_list clone() const {
				node* currentNode = node0;
				node* newNode0 = new node(currentNode->value);
				node* newNodeF = newNode0;
				for (int i = 1; i < length; i++) {
					currentNode = currentNode->next;
					newNodeF = (newNode0->next = new node(currentNode->value));
				}
				return linked_list(newNode0, newNodeF, length);
			}

			void append(T item) {
				nodeF = (nodeF->next = new node(item));
				length += 1;
			}

			void append(T* items, int copyLength) {
				for (int i = 0; i < copyLength; i++) {
					nodeF = (nodeF->next = new node(items[i]));
				}
				length += copyLength;
			}

			void append(const linked_list& items) {
				node* currentNode = items.node0;
				for (int i = 0; i < items.length; i++) {
					nodeF = (nodeF->next = new node(currentNode->value));
					currentNode = currentNode->next;
				}
				length += items.length;
			}

			void operator +=(T item) {
				append(item);
			}

			void operator +=(const linked_list& items) {
				append(items);
			}

			void insertA(T item, int i) { //above
				node* lowerNode;
				if (i) {
					node* upperNode = (*node0)[i - 1];
					lowerNode = upperNode->next;
					upperNode->next = new node(item);
					upperNode->next->next = lowerNode;
				}
				else {
					lowerNode = node0;
					node0 = new node(item);
					node0->next = lowerNode;
				}
				length += 1;
			}

			void insertB(T item, int i) { //below
				if (length == (i + 1)) {
					nodeF = (nodeF->next = new node(item));
				}
				else {
					node* upperNode;
					node* lowerNode;
					upperNode = (*node0)[i];
					lowerNode = upperNode->next;
					upperNode->next = new node(item);
					upperNode->next->next = lowerNode;
				}
				length += 1;
			}

			void insertB(T item, void* pos) {
				node* upperNode = reinterpret_cast<node*>(pos);
				if (upperNode == nodeF) {
					nodeF = (nodeF->next = new node(item));
				}
				else {
					node* lowerNode = upperNode->next;
					upperNode->next = new node(item);
					upperNode->next->next = lowerNode;
				}
				length += 1;
			}

			void remove(int i) {
				if (i) {
					node* upperNode = (*node0)[i - 1];
					node* removedNode = upperNode->next;
					node* lowerNode = removedNode->next;
					removedNode->next = 0;
					delete removedNode;
					upperNode->next = lowerNode;
					length -= 1;
					if (length == i) {
						nodeF = upperNode;
					}
				}
				else {
					if (length > 1) {
						node* removedNode = node0;
						node* lowerNode = removedNode->next;
						removedNode->next = 0;
						delete removedNode;
						node0 = lowerNode;
						length -= 1;
						if (length == i) {
							nodeF = lowerNode;
						}
					}
					else {
						throw 0;
					}
				}
			}

			void remove(void* nodePosPtr) {
				node* nodePos = reinterpret_cast<node*>(nodePosPtr);
				if (nodePos == node0) {
					if (length > 1) {
						node* removedNode = node0;
						node* lowerNode = removedNode->next;
						removedNode->next = 0;
						delete removedNode;
						node0 = lowerNode;
						length -= 1;
						if (length == 1) {
							nodeF = node0;
						}
					}
					else {
						throw 0;
					}
				}
				else {
					node* upperNode = node0;
					node* removedNode;
					node* lowerNode;
					for (int i = 0; i < length; i++) {
						if (nodePos == upperNode->next) {
							removedNode = upperNode->next;
							lowerNode = removedNode->next;
							removedNode->next = 0;
							delete removedNode;
							upperNode->next = lowerNode;
							length -= 1;
							if (length == (i + 1)) {
								nodeF = upperNode;
							}
							break;
						}
						else {
							upperNode = upperNode->next;
						}
					}
				}
			}

			T& operator [](int i) {
				if (i >= length) { throw 0; }
				return (*node0)[i]->value;
			}

			T& operator [](void* pos) {
				return reinterpret_cast<node*>(pos)->value;
			}

			int size() {
				return length;
			}

			~linked_list() {
				delete_node(node0);
			}
	};

	template <typename T>
	std::ostream& operator << (std::ostream& os, linked_list<T>& list) {
		os << "{ ";
		int lastIndex = list.size() - 1;
		void* node = list.getStart();
		for (int i = 0; i < lastIndex; i++) {
			os << list[node] << ", ";
			node = list.next(node);
		}
		os << list[node] << " }";
		return os;
	}

	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////

	template <typename K, typename V>
	struct binary_search_tree {
		V value;
		K key;
		int size;
		binary_search_tree* branches[2];
		binary_search_tree(K newKey, V newValue) {
			key = newKey;
			value = newValue;
			branches[0] = 0;
			branches[1] = 0;
			size = 1;
		}
		binary_search_tree(const binary_search_tree& bst) {
			value = bst.value;
			key = bst.key;
			size = bst.size;
			branches[0] = bst.branches[0];
			branches[1] = bst.branches[1];
		}
		V& safeSearch(K searchKey, bool* success) {
			binary_search_tree* currentBranch = this;
			while (currentBranch) {
				K currentKey = currentBranch->key;
				if (searchKey == currentKey) {
					*success = true;
					return currentBranch->value;
				}
				else {
					currentBranch = currentBranch->branches[searchKey > currentKey];
				}
			}
			*success = false;
			V defaultValue;
			return defaultValue;
		}
		V& operator [](K searchKey) {
			binary_search_tree* currentBranch = this;
			while (currentBranch) {
				K currentKey = currentBranch->key;
				if (searchKey == currentKey) {
					return currentBranch->value;
				}
				else {
					currentBranch = currentBranch->branches[searchKey > currentKey];
				}
			}
			throw 0;
		}
		~binary_search_tree() {
			delete branches[0];
			delete branches[1];
		}
	};

	template <typename K, typename V>
	V search(binary_search_tree<K, V>* bst, K key) {
		while (!(key == bst->key)) {
			bst = bst->branches[key > bst->key];
		}
		return bst->value;
	}

	template <typename K, typename V>
	V ssearch(binary_search_tree<K, V>* bst, K key) {
		while (true) {
			if (bst) {
				if (!(key == bst->key)) {
					bst = bst->branches[key > bst->key];
				}
				else { return bst->value; }
			}
			else { throw 0; }
		}
	}

	template <typename K, typename V>
	void insert(binary_search_tree<K, V>* bst, K key, V value, bool* extended) {
		// recursive version
		if (key == bst->key) { bst->value = value; }
		else {
			bool nextBranchIndex = key > bst->key;
			binary_search_tree<K, V>* nextBranch = bst->branches[nextBranchIndex];
			if (nextBranch) {
				insert(nextBranch, key, value);
			}
			else {
				bst->branches[nextBranchIndex] = new binary_search_tree(key, value);
				*extended = true;
			}
		}
		if (*extended) { bst->size += 1; }
	}

	template <typename K, typename V>
	void insertnr(binary_search_tree<K, V>* bst, K key, V value) {
		// non-recursive version
		bool extended = false;
		binary_search_tree<K, V>* currentBranch = bst;
		binary_search_tree<K, V>* nextBranch;
		bool nextBranchIndex;
		int maxLevels = bst->size;
		while (maxLevels--) {
			if (key == currentBranch->key) {
				currentBranch->value = value;
				maxLevels = 0;
				break;
			}
			else {
				nextBranchIndex = key > currentBranch->key;
				nextBranch = currentBranch->branches[nextBranchIndex];
				if (nextBranch) {
					currentBranch = nextBranch;
				}
				else {
					currentBranch->branches[nextBranchIndex] = new binary_search_tree(key, value);
					extended = true;
					maxLevels = 0;
					break;
				}
			}
		}
		if (extended) {
			currentBranch = bst;
			maxLevels = bst->size;
			while (maxLevels--) {
				if (key == currentBranch->key) { maxLevels = 0; break; }
				else {
					nextBranchIndex = key > currentBranch->key;
					nextBranch = currentBranch->branches[nextBranchIndex];
					currentBranch->size += 1;
					if (nextBranch) {
						currentBranch = nextBranch;
					}
					else {
						maxLevels = 0;
						break;
					}
				}
			}
		}
	}

	template <typename K, typename V>
	void insert(binary_search_tree<K, V>* bst, K key, V value) {
		bool extended = false;
		//insert(bst, key, value, &extended);
		insertnr(bst, key, value);
	}

	template <typename K, typename V>
	void insertExtreme(binary_search_tree<K, V>* bst, binary_search_tree<K, V>* extremeBST, bool side) {
		int extremeSize = extremeBST->size;
		while (true) {
			bst->size += extremeSize;
			if (bst->branches[side]) { bst = bst->branches[side]; }
			else { bst->branches[side] = extremeBST; break; }
		}
	}

	template <typename K, typename V>
	void rotate(binary_search_tree<K, V>* bst, bool rotationIndex) {
		//form branch from top to be inserted as an extreme and bumped down into the tree
		binary_search_tree<K, V>* bumpedDown = new binary_search_tree<K, V>(bst->key, bst->value);
		bumpedDown->branches[!rotationIndex] = bst->branches[!rotationIndex];
		if (bumpedDown->branches[!rotationIndex]) {
			bumpedDown->size += bumpedDown->branches[!rotationIndex]->size;
		}

		//copy rotated node to top of tree
		binary_search_tree<K, V>* rotatedNode = bst->branches[rotationIndex];
		bst->key = rotatedNode->key;
		bst->value = rotatedNode->value;
		bst->size = rotatedNode->size;
		bst->branches[1] = rotatedNode->branches[1];
		bst->branches[0] = rotatedNode->branches[0];

		//delete out of use rotated node
		rotatedNode->branches[0] = 0;
		rotatedNode->branches[1] = 0;
		delete rotatedNode;

		//bump down the extreme branch
		int extremeSize = bumpedDown->size;
		bool extremeSide = !rotationIndex;
		binary_search_tree<K, V>* currentBranch = bst;
		while (true) {
			currentBranch->size += extremeSize;
			if (currentBranch->branches[extremeSide]) {
				currentBranch = currentBranch->branches[extremeSide];
			}
			else {
				currentBranch->branches[extremeSide] = bumpedDown;
				break;
			}
		}
	}

	template <typename K, typename V>
	void balance(binary_search_tree<K, V>* bst, int sensitivity) {
		//std::cout << "STARTING BRANCH LEVEL" << std::endl;
		bool b0exists;
		bool b1exists;
		int b0size;
		int b1size;
		int difference;
		bool rotation;
		while (true) {
			b0exists = bst->branches[0];
			b1exists = bst->branches[1];
			if (b0exists) { b0size = bst->branches[0]->size; }
			else { b0size = 0; }
			if (b1exists) { b1size = bst->branches[1]->size; }
			else { b1size = 0; }
			difference = b1size - b0size;
			//std::cout << "diff " << difference << std::endl;
			//std::cout << "sizes " << b0size << ", " << b1size << std::endl;
			if (std::abs(difference) > (1 + sensitivity)) {
				rotation = difference > 0;
				rotate(bst, rotation);
			}
			else { break; }
		}
		if (b0exists) { balance(bst->branches[0], sensitivity); }
		if (b1exists) { balance(bst->branches[1], sensitivity); }
	}

	template <typename K, typename V>
	void remove(binary_search_tree<K, V>* bst, K key) {
		//locating the node to be removed with safe search (throws int)
		while (true) {
			if (bst) {
				bst->size -= 1;
				if (key == bst->key) { break; }
				else { bst = bst->branches[key > bst->key]; }
			}
			else { throw 0; }
		}

		//evaluating the existance and size of branches
		bool b0exists = bst->branches[0];
		bool b1exists = bst->branches[1];
		int b0size, b1size;
		if (b0exists) { b0size = bst->branches[0]->size; }
		else { b0size = 0; }
		if (b1exists) { b1size = bst->branches[1]->size; }
		else { b1size = 0; }

		if (b0size || b1size) {
			//compacting the contents of lower branches into one larger branch that will replace the removed node
			bool greaterBranchIdx = b1size > b0size;
			binary_search_tree<K, V>* greaterBranch = bst->branches[greaterBranchIdx];
			if (bst->branches[!greaterBranchIdx]) {
				insertExtreme(greaterBranch, bst->branches[!greaterBranchIdx], !greaterBranchIdx);
			}

			//copying larger branch's top node's data to the removed node's location
			bst->key = greaterBranch->key;
			bst->value = greaterBranch->value;
			bst->branches[0] = greaterBranch->branches[0];
			bst->branches[1] = greaterBranch->branches[1];

			//delete the node at the top of the larger branch
			greaterBranch->branches[0] = 0;
			greaterBranch->branches[1] = 0;
			delete greaterBranch;
		}
		else {
			delete bst;
		}
	}

	template <typename K, typename V>
	void bremove(binary_search_tree<K, V>* bst, K key, int sensitivity) {
		//locating the node to be removed with safe search (throws int)
		while (true) {
			if (bst) {
				bst->size -= 1;
				if (key == bst->key) { break; }
				else { bst = bst->branches[key > bst->key]; }
			}
			else { throw 0; }
		}

		//evaluating the existance and size of branches
		bool b0exists = bst->branches[0];
		bool b1exists = bst->branches[1];
		int b0size, b1size;
		if (b0exists) { b0size = bst->branches[0]->size; }
		else { b0size = 0; }
		if (b1exists) { b1size = bst->branches[1]->size; }
		else { b1size = 0; }

		if (b0size || b1size) {
			//compacting the contents of lower branches into one larger branch that will replace the removed node
			bool greaterBranchIdx = b1size > b0size;
			binary_search_tree<K, V>* greaterBranch = bst->branches[greaterBranchIdx];
			if (bst->branches[!greaterBranchIdx]) {
				insertExtreme(greaterBranch, bst->branches[!greaterBranchIdx], !greaterBranchIdx);
				balance(greaterBranch, sensitivity);
			}

			//copying larger branch's top node's data to the removed node's location
			bst->key = greaterBranch->key;
			bst->value = greaterBranch->value;
			bst->branches[0] = greaterBranch->branches[0];
			bst->branches[1] = greaterBranch->branches[1];

			//delete the node at the top of the larger branch
			greaterBranch->branches[0] = 0;
			greaterBranch->branches[1] = 0;
			delete greaterBranch;
		}
		else {
			delete bst;
		}
	}

	template <typename K, typename V>
	struct ordered_array {
		V* values;
		K* keys;
		int len;
		ordered_array(const ordered_array& oarr) {
			values = new V[len = oarr.len];
			keys = new K[len = oarr.len];
			for (int i = 0; i < len; i++) {
				values[i] = oarr.values[i];
				keys[i] = oarr.keys[i];
			}
		}
		void add(const binary_search_tree<K, V>* bst, int* idx) {
			if (bst->branches[0]) add(bst->branches[0], idx);
			values[*idx] = bst->value;
			keys[(*idx)++] = bst->key;
			if (bst->branches[1]) add(bst->branches[1], idx);
		}
		ordered_array(const binary_search_tree<K, V>* bst) {
			values = new V[len = bst->size];
			keys = new K[len = bst->size];
			int idx = 0;
			add(bst, &idx);
		}
		ordered_array operator =(const binary_search_tree <K, V>* bst) {
			values = new V[len = bst->size];
			keys = new K[len = bst->size];
			int idx = 0;
			add(bst, &idx);
			return *this;
		}
		V& operator [](int i) {
			return values[i];
		}
		~ordered_array() {
			delete[] values;
			delete[] keys;
		}
	};

	template<typename K, typename V>
	int isearch(const ordered_array<K, V>* oarr, K key) {
		int marker = -1;
		int direction = 1;
		bool continuous = true;
		unsigned int range;
		for (unsigned int size = oarr->len; size > 0; size = range - continuous * (1 + 2 * range - size)) {
			range = size >> 1;
			marker += direction * (range + 1);
			if (key == oarr->keys[marker]) return marker;
			continuous = direction == (2 * (key > oarr->keys[marker]) - 1);
			direction *= 2 * continuous - 1;
		}
		return marker;
	}

	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////

	// linkarray_list, array_list, queue (rotating array_list), matrix
}