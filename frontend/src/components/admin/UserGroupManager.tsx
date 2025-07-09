import React, { useState, useEffect } from 'react';
import { FaUsers, FaUserCog, FaTrash, FaEdit, FaSave, FaTimes, FaHammer, FaCheck, FaTimes as FaX } from 'react-icons/fa';

interface User {
  user_id: string;
  email: string;
  group_name: string;
  group_description?: string;
  created_at: string;
}

interface Group {
  group_id: string;
  name: string;
  description?: string;
  created_at: string;
}

interface Workflow {
  workflow_id: number;
  title: string;
  description?: string;
  initial_prompt?: string;
}

interface WorkflowAccess {
  group_id: string;
  workflow_id: number;
  has_access: boolean;
}

const UserGroupManager: React.FC = () => {
  const [users, setUsers] = useState<User[]>([]);
  const [groups, setGroups] = useState<Group[]>([]);
  const [workflows, setWorkflows] = useState<Workflow[]>([]);
  const [workflowAccess, setWorkflowAccess] = useState<WorkflowAccess[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [editingGroup, setEditingGroup] = useState<string | null>(null);
  const [newGroup, setNewGroup] = useState({ name: '', description: '' });
  const [showCreateForm, setShowCreateForm] = useState(false);
  const [activeTab, setActiveTab] = useState<'users' | 'workflows'>('users');

  // Function to sort groups in desired order: admin, researcher, basic
  const sortGroups = (groups: Group[]) => {
    const order = ['admin', 'researcher', 'basic'];
    return groups.sort((a, b) => {
      const aIndex = order.indexOf(a.name);
      const bIndex = order.indexOf(b.name);
      // If both groups are in the predefined order, sort by their position
      if (aIndex !== -1 && bIndex !== -1) {
        return aIndex - bIndex;
      }
      // If only one is in the predefined order, prioritize it
      if (aIndex !== -1) return -1;
      if (bIndex !== -1) return 1;
      // If neither is in the predefined order, sort alphabetically
      return a.name.localeCompare(b.name);
    });
  };

  // Function to sort workflows according to default_workflows.json order
  const sortWorkflows = (workflows: Workflow[]) => {
    const defaultOrder = [1, 2, 3, 4, 5]; // workflow_ids in order from default_workflows.json
    return workflows.sort((a, b) => {
      const aIndex = defaultOrder.indexOf(a.workflow_id);
      const bIndex = defaultOrder.indexOf(b.workflow_id);
      // If both workflows are in the predefined order, sort by their position
      if (aIndex !== -1 && bIndex !== -1) {
        return aIndex - bIndex;
      }
      // If only one is in the predefined order, prioritize it
      if (aIndex !== -1) return -1;
      if (bIndex !== -1) return 1;
      // If neither is in the predefined order, sort by workflow_id
      return a.workflow_id - b.workflow_id;
    });
  };

  // Load data on component mount
  useEffect(() => {
    loadData();
  }, []);

  const loadData = async () => {
    try {
      setLoading(true);
      const [usersResponse, groupsResponse] = await Promise.all([
        fetch('/api/admin/users'),
        fetch('/api/admin/groups')
      ]);

      if (!usersResponse.ok || !groupsResponse.ok) {
        throw new Error('Failed to load data');
      }

      const usersData = await usersResponse.json();
      const groupsData = await groupsResponse.json();

      setUsers(usersData.users);
      setGroups(groupsData.groups);

      // Load workflow access data (this includes all workflows)
      await loadWorkflowAccess();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An error occurred');
    } finally {
      setLoading(false);
    }
  };

  const loadWorkflowAccess = async () => {
    try {
      const response = await fetch('/api/admin/workflow-access');
      if (!response.ok) {
        throw new Error('Failed to load workflow access data');
      }
      
      const data = await response.json();
      console.log('Workflow access data:', data);
      setWorkflowAccess(data.access_matrix);
      // Set workflows from the admin endpoint (includes all workflows)
      setWorkflows(data.workflows || []);
    } catch (err) {
      console.error('Failed to load workflow access:', err);
    }
  };

  const toggleWorkflowAccess = async (groupId: string, workflowId: number, currentAccess: boolean) => {
    try {
      const method = currentAccess ? 'DELETE' : 'POST';
      const response = await fetch(`/api/admin/groups/${groupId}/workflows/${workflowId}/access`, {
        method,
        headers: { 'Content-Type': 'application/json' }
      });

      if (!response.ok) {
        throw new Error(`Failed to ${currentAccess ? 'revoke' : 'grant'} workflow access`);
      }

      // Update local state
      setWorkflowAccess(prev => 
        prev.map(access => 
          access.group_id === groupId && access.workflow_id === workflowId
            ? { ...access, has_access: !currentAccess }
            : access
        )
      );
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update workflow access');
    }
  };

  const updateUserGroup = async (userId: string, groupId: string) => {
    try {
      const response = await fetch(`/api/admin/users/${userId}/group`, {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ group_id: groupId })
      });

      if (!response.ok) {
        throw new Error('Failed to update user group');
      }

      // Reload data to reflect changes
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update user group');
    }
  };

  const createGroup = async () => {
    try {
      const response = await fetch('/api/admin/groups', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(newGroup)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || 'Failed to create group');
      }

      setNewGroup({ name: '', description: '' });
      setShowCreateForm(false);
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create group');
    }
  };

  const updateGroup = async (groupId: string, name: string, description: string) => {
    try {
      const response = await fetch(`/api/admin/groups/${groupId}`, {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name, description })
      });

      if (!response.ok) {
        throw new Error('Failed to update group');
      }

      setEditingGroup(null);
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update group');
    }
  };

  const deleteGroup = async (groupId: string, groupName: string) => {
    if (!confirm(`Are you sure you want to delete the group "${groupName}"? Users in this group will be moved to the basic group.`)) {
      return;
    }

    try {
      const response = await fetch(`/api/admin/groups/${groupId}`, {
        method: 'DELETE'
      });

      if (!response.ok) {
        throw new Error('Failed to delete group');
      }

      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to delete group');
    }
  };

  if (loading) {
    return (
      <div className="flex justify-center items-center h-64">
        <div className="text-white">Loading...</div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="bg-red-600 text-white p-4 rounded-lg mb-4">
        <div className="font-bold">Error:</div>
        <div>{error}</div>
        <button 
          onClick={() => setError(null)}
          className="mt-2 px-3 py-1 bg-red-700 rounded hover:bg-red-800"
        >
          Dismiss
        </button>
      </div>
    );
  }

  return (
    <div className="max-w-6xl mx-auto p-6">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-white mb-2 flex items-center gap-3">
          <FaUsers className="text-blue-400" />
          User Group Management
        </h1>
        <p className="text-gray-300">Manage user groups, assign users to different permission levels, and control workflow access permissions for each group. Create custom groups, modify existing ones, and ensure proper access control across the platform.</p>
      </div>

      {/* Tab Navigation */}
      <div className="flex border-b border-gray-700 mb-6">
        <button
          onClick={() => setActiveTab('users')}
          className={`px-6 py-3 font-medium rounded-t-lg transition-colors ${
            activeTab === 'users'
              ? 'bg-gray-800 text-white border-b-2 border-blue-400'
              : 'text-gray-400 hover:text-white hover:bg-gray-700'
          }`}
        >
          <FaUsers className="inline mr-2" />
          User Management
        </button>
        <button
          onClick={() => setActiveTab('workflows')}
          className={`px-6 py-3 font-medium rounded-t-lg transition-colors ${
            activeTab === 'workflows'
              ? 'bg-gray-800 text-white border-b-2 border-blue-400'
              : 'text-gray-400 hover:text-white hover:bg-gray-700'
          }`}
        >
          <FaHammer className="inline mr-2" />
          Workflow Access
        </button>
      </div>

      {/* User Management Tab */}
      {activeTab === 'users' && (
        <>
          {/* Groups Section */}
          <div className="bg-gray-800 rounded-lg p-6 mb-8">
            <div className="flex justify-between items-center mb-6">
              <h2 className="text-xl font-semibold text-white flex items-center gap-2">
                <FaUserCog className="text-blue-400" />
                User Groups
              </h2>
              <button
                onClick={() => setShowCreateForm(true)}
                className="bg-blue-600 hover:bg-blue-700 text-white px-4 py-2 rounded-lg flex items-center gap-2"
              >
                Create Group
              </button>
            </div>

            {/* Create Group Form */}
            {showCreateForm && (
              <div className="bg-gray-700 rounded-lg p-4 mb-4">
                <h3 className="text-lg font-medium text-white mb-3">Create New Group</h3>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
                  <div>
                    <label className="block text-sm font-medium text-gray-300 mb-1">Name</label>
                    <input
                      type="text"
                      value={newGroup.name}
                      onChange={(e) => setNewGroup({ ...newGroup, name: e.target.value })}
                      className="w-full px-3 py-2 bg-gray-600 border border-gray-500 rounded-md text-white focus:outline-none focus:border-blue-500"
                      placeholder="Enter group name"
                    />
                  </div>
                  <div>
                    <label className="block text-sm font-medium text-gray-300 mb-1">Description</label>
                    <input
                      type="text"
                      value={newGroup.description}
                      onChange={(e) => setNewGroup({ ...newGroup, description: e.target.value })}
                      className="w-full px-3 py-2 bg-gray-600 border border-gray-500 rounded-md text-white focus:outline-none focus:border-blue-500"
                      placeholder="Enter description"
                    />
                  </div>
                </div>
                <div className="flex gap-2">
                  <button
                    onClick={createGroup}
                    disabled={!newGroup.name.trim()}
                    className="bg-green-600 hover:bg-green-700 disabled:bg-gray-600 text-white px-4 py-2 rounded-lg flex items-center gap-2"
                  >
                    <FaSave />
                    Create
                  </button>
                  <button
                    onClick={() => {
                      setShowCreateForm(false);
                      setNewGroup({ name: '', description: '' });
                    }}
                    className="bg-gray-600 hover:bg-gray-700 text-white px-4 py-2 rounded-lg flex items-center gap-2"
                  >
                    <FaTimes />
                    Cancel
                  </button>
                </div>
              </div>
            )}

            {/* Groups List */}
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
              {groups.map((group) => (
                <div key={group.group_id} className="bg-gray-700 rounded-lg p-4">
                  {editingGroup === group.group_id ? (
                    <GroupEditForm
                      group={group}
                      onSave={(name, description) => updateGroup(group.group_id, name, description)}
                      onCancel={() => setEditingGroup(null)}
                    />
                  ) : (
                    <div>
                      <div className="flex justify-between items-start mb-2">
                        <h3 className="text-lg font-medium text-white">{group.name}</h3>
                        <div className="flex gap-1">
                          <button
                            onClick={() => setEditingGroup(group.group_id)}
                            className="text-blue-400 hover:text-blue-300 p-1"
                            title="Edit group"
                          >
                            <FaEdit />
                          </button>
                          {group.name !== 'admin' && (
                            <button
                              onClick={() => deleteGroup(group.group_id, group.name)}
                              className="text-red-400 hover:text-red-300 p-1"
                              title="Delete group"
                            >
                              <FaTrash />
                            </button>
                          )}
                        </div>
                      </div>
                      <p className="text-gray-300 text-sm mb-3">
                        {group.description || 'No description'}
                      </p>
                      <div className="text-xs text-gray-400">
                        Users: {users.filter(u => u.group_name === group.name).length}
                      </div>
                    </div>
                  )}
                </div>
              ))}
            </div>
          </div>

          {/* Users Section */}
          <div className="bg-gray-800 rounded-lg p-6">
            <h2 className="text-xl font-semibold text-white mb-6 flex items-center gap-2">
              <FaUsers className="text-blue-400" />
              Users
            </h2>
            
            <div className="overflow-x-auto">
              <table className="w-full text-left">
                <thead>
                  <tr className="border-b border-gray-700">
                    <th className="py-3 px-4 text-gray-300 font-medium">Email</th>
                    <th className="py-3 px-4 text-gray-300 font-medium">Group</th>
                    <th className="py-3 px-4 text-gray-300 font-medium">Created</th>
                    <th className="py-3 px-4 text-gray-300 font-medium">Actions</th>
                  </tr>
                </thead>
                <tbody>
                  {users.map((user) => (
                    <tr key={user.user_id} className="border-b border-gray-700 hover:bg-gray-700">
                      <td className="py-3 px-4 text-white">{user.email}</td>
                      <td className="py-3 px-4">
                        <select
                          value={groups.find(g => g.name === user.group_name)?.group_id || ''}
                          onChange={(e) => updateUserGroup(user.user_id, e.target.value)}
                          className="bg-gray-600 border border-gray-500 rounded px-2 py-1 text-white text-sm focus:outline-none focus:border-blue-500"
                        >
                          {groups.map((group) => (
                            <option key={group.group_id} value={group.group_id}>
                              {group.name}
                            </option>
                          ))}
                        </select>
                      </td>
                      <td className="py-3 px-4 text-gray-300 text-sm">
                        {new Date(user.created_at).toLocaleDateString()}
                      </td>
                      <td className="py-3 px-4">
                        <span className="text-xs text-gray-400">
                          {user.group_description || 'No description'}
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </>
      )}

      {/* Workflow Access Management Tab */}
      {activeTab === 'workflows' && (
        <div className="bg-gray-800 rounded-lg p-6">
          <div className="mb-6">
            <h2 className="text-xl font-semibold text-white flex items-center gap-2">
              <FaHammer className="text-blue-400" />
              Workflow Access Management
            </h2>
            <p className="text-gray-300 mt-2">Control which user groups can access which workflows</p>
          </div>

          {workflows.length > 0 && groups.length > 0 ? (
            <div className="overflow-x-auto">
              <table className="w-full text-left">
                <thead>
                  <tr className="border-b border-gray-700">
                    <th className="py-3 px-4 text-gray-300 font-medium">Workflow</th>
                    {sortGroups(groups).map((group) => (
                      <th key={group.group_id} className="py-3 px-4 text-gray-300 font-medium text-center">
                        {group.name}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {sortWorkflows(workflows).map((workflow) => (
                    <tr key={workflow.workflow_id} className="border-b border-gray-700 hover:bg-gray-700">
                      <td className="py-3 px-4">
                        <div className="text-white font-medium">{workflow.title}</div>
                        <div className="text-gray-400 text-sm">{workflow.description}</div>
                      </td>
                      {sortGroups(groups).map((group) => {
                        const access = workflowAccess.find(
                          a => a.group_id === group.group_id && a.workflow_id === workflow.workflow_id
                        );
                        const hasAccess = access?.has_access || false;
                        
                        console.log(`Group ${group.name} (${group.group_id}) - Workflow ${workflow.title} (${workflow.workflow_id}): hasAccess=${hasAccess}`, access);
                        console.log('All workflow access data:', workflowAccess);
                        
                        return (
                          <td key={group.group_id} className="py-3 px-4 text-center">
                            <button
                              onClick={() => toggleWorkflowAccess(group.group_id, workflow.workflow_id, hasAccess)}
                              className={`w-12 h-12 rounded-full flex items-center justify-center transition-colors border-2 ${
                                hasAccess
                                  ? 'bg-green-600 hover:bg-green-700 text-white border-green-500'
                                  : 'bg-gray-700 hover:bg-gray-600 text-gray-300 border-gray-500'
                              }`}
                              title={`${hasAccess ? 'Revoke' : 'Grant'} access for ${group.name}`}
                            >
                              <span className="text-3xl font-bold">
                                {hasAccess ? '✓' : '✗'}
                              </span>
                            </button>
                          </td>
                        );
                      })}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          ) : (
            <div className="text-center py-8">
              <div className="text-gray-400">No workflows or groups available</div>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

// Group Edit Form Component
interface GroupEditFormProps {
  group: Group;
  onSave: (name: string, description: string) => void;
  onCancel: () => void;
}

const GroupEditForm: React.FC<GroupEditFormProps> = ({ group, onSave, onCancel }) => {
  const [name, setName] = useState(group.name);
  const [description, setDescription] = useState(group.description || '');

  return (
    <div>
      <div className="mb-3">
        <label className="block text-sm font-medium text-gray-300 mb-1">Name</label>
        <input
          type="text"
          value={name}
          onChange={(e) => setName(e.target.value)}
          className="w-full px-3 py-2 bg-gray-600 border border-gray-500 rounded-md text-white focus:outline-none focus:border-blue-500"
        />
      </div>
      <div className="mb-3">
        <label className="block text-sm font-medium text-gray-300 mb-1">Description</label>
        <input
          type="text"
          value={description}
          onChange={(e) => setDescription(e.target.value)}
          className="w-full px-3 py-2 bg-gray-600 border border-gray-500 rounded-md text-white focus:outline-none focus:border-blue-500"
        />
      </div>
      <div className="flex gap-2">
        <button
          onClick={() => onSave(name, description)}
          disabled={!name.trim()}
          className="bg-green-600 hover:bg-green-700 disabled:bg-gray-600 text-white px-3 py-1 rounded text-sm flex items-center gap-1"
        >
          <FaSave />
          Save
        </button>
        <button
          onClick={onCancel}
          className="bg-gray-600 hover:bg-gray-700 text-white px-3 py-1 rounded text-sm flex items-center gap-1"
        >
          <FaTimes />
          Cancel
        </button>
      </div>
    </div>
  );
};

export default UserGroupManager; 