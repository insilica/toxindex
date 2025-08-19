import React, { createContext, useContext, useState, useCallback, useEffect } from "react";

export type ChatSessionContextType = {
  chatSessions: any[];
  selectedSessionId: string | null;
  setSelectedSessionId: (id: string | null) => void;
  refetchChatSessions: () => Promise<void>;
};

const ChatSessionContext = createContext<ChatSessionContextType | undefined>(undefined);

export const ChatSessionProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  console.log('[ChatSessionProvider] ChatSessionProvider mounting');
  const [chatSessions, setChatSessions] = useState<any[]>([]);
  const [selectedSessionId, setSelectedSessionId] = useState<string | null>(null);

  useEffect(() => {
    console.log('[ChatSessionProvider] ChatSessionProvider mounted');
    return () => {
      console.log('[ChatSessionProvider] ChatSessionProvider unmounting');
    };
  }, []);

  const refetchChatSessions = useCallback(async () => {
    console.log('[ChatSessionProvider] Refetching chat sessions...');
    const res = await fetch(`/api/chat_sessions`, { credentials: 'include' });
    const data = await res.json();
    console.log('[ChatSessionProvider] Chat sessions fetched:', data.sessions);
    setChatSessions(data.sessions || []);
    if (data.sessions && data.sessions.length > 0) {
      console.log('[ChatSessionProvider] Setting selected session to first session:', data.sessions[0].session_id);
      setSelectedSessionId(data.sessions[0].session_id);
    } else {
      console.log('[ChatSessionProvider] No sessions found, clearing selected session');
      setSelectedSessionId(null);
    }
  }, []);

  console.log('[ChatSessionProvider] ChatSessionProvider rendering with state:', {
    chatSessionsCount: chatSessions.length,
    selectedSessionId
  });

  return (
    <ChatSessionContext.Provider value={{
      chatSessions,
      selectedSessionId,
      setSelectedSessionId,
      refetchChatSessions
    }}>
      {children}
    </ChatSessionContext.Provider>
  );
};

export const useChatSession = () => {
  const context = useContext(ChatSessionContext);
  if (!context) throw new Error("useChatSession must be used within a ChatSessionProvider");
  return context;
}; 